#!/usr/bin/env python3
"""
tree_matrix.py - Gene x species state matrix, aligned on the phylogeny.

This is the DISPLAY layer of the chain, and the only place where every state is
shown as it was assigned. Nothing is regrouped, nothing is thresholded, nothing
is decided:

    build_candidates.py   assigns one state per gene x species (states.tsv)
    tree_matrix.py        shows those states, all of them, side by side  <- here
    build_ordered_matrix  regroups them into presence / absence / unscorable
    select_nested.py      decides

That separation is the point. 'coding_flank' in particular is visible here as
itself, distinct from 'absent': a neighbour with coding sequence past the
alignment is not a neighbour without a homolog, and the two must not look alike
at the stage where one is still looking rather than concluding.

States drawn, as produced by build_candidates.py:

  absent            no homolog passed the filters. Says nothing
  coding_flank      the homolog has coding sequence past the alignment, so it
                    cannot testify about the non-coding region. Undecidable, and
                    undecidable for the whole species once aggregation puts
                    coding_flank first (--coding-flank-priority high)
  complete          the species already carries the extension in an annotated
                    protein. Argues against a recent elongation
  no_trace          small category, nothing recovered in the flank
  not_analysable    N-ter only, and only with --separate-not-analysable: the
                    upstream region could not be read at all
  annotation_doubt  N-ter only: a usable in-frame ATG sits upstream with no stop
                    before the annotated start, so the neighbour could translate
                    that region and may simply be annotated from the wrong ATG.
                    The C-ter cannot have this state, a stop being unique per
                    frame
  weak_trace        small category, ground gained but not a single contiguous
                    alignment over most of the extension
  elongation        small category, the trace is there: one local alignment
                    covering most of the extension and reaching the conserved
                    core

The focal species is drawn as its own row and its own colour: it carries the
extension by construction, so it has no state in states.tsv and must not be
painted as 'absent'.

No filtering is applied unless asked. Every gene in states.tsv is drawn, or the
subset given by --genes-from, which takes the list select_nested.py writes.

Usage:
    python3 tree_matrix.py --side nter
    python3 tree_matrix.py --side nter --order nested
    python3 tree_matrix.py --side nter --genes-from top_nter.txt
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import polars as pl
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
from Bio import Phylo


# Per side: state name -> (integer code, colour, legend label).
# Codes are display identities only. The conservative aggregation of several
# subjects of one species into a single state happens in build_candidates.py and
# nowhere else; duplicating it here is how two definitions drift apart.
STATES_CTER = {
    "absent":            (0, "#f2f2f2", "aucun homologue"),
    "coding_flank":      (1, "#9ba7b0", "séquence codante en aval, indécidable"),
    "complete":          (2, "#2f4b7c", "extension déjà présente"),
    "no_trace":          (3, "#f6c85f", "aucune trace en aval"),
    "weak_trace":        (4, "#e8a33d", "trace partielle"),
    "elongation":        (5, "#c3423f", "trace retrouvée"),
    "focal":             (6, "#ffffff", "espèce focale"),
}

STATES_NTER = {
    "absent":            (0, "#f2f2f2", "aucun homologue"),
    "coding_flank":      (1, "#9ba7b0", "séquence codante en amont, indécidable"),
    "complete":          (2, "#2f4b7c", "extension déjà présente"),
    "no_trace":          (3, "#f6c85f", "aucune trace en amont"),
    "not_analysable":    (4, "#cdbfe0", "région amont illisible"),
    "annotation_doubt":  (5, "#6f9bd1", "début d'annotation douteux"),
    "weak_trace":        (6, "#e8a33d", "trace partielle"),
    "elongation":        (7, "#c3423f", "trace retrouvée"),
    "focal":             (8, "#ffffff", "espèce focale"),
}


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--side", choices=["cter", "nter"], default="nter")
    p.add_argument("--states", type=Path,
                   default=Path("analysis/tables/states.tsv"),
                   help="states table from build_candidates.py")
    p.add_argument("--pairs", type=Path,
                   default=Path("analysis/tables/pairs.tsv"),
                   help="pairs table, used only to label columns by locus")
    p.add_argument("--nested", type=Path, default=None,
                   help="nested table from select_nested.py (default: "
                        "analysis/tables/{side}_nested.tsv). Only read by "
                        "--order nested")
    p.add_argument("--genes-from", type=Path, default=None,
                   help="file with one qseqid per line, as written by "
                        "select_nested.py --out-list; restricts the columns")
    p.add_argument("--cladogram", action="store_true", default=True,
                   help="space nodes evenly and align the tips (default)")
    p.add_argument("--true-lengths", dest="cladogram", action="store_false",
                   help="respect branch lengths instead")
    p.add_argument("--outdir", type=Path, default=Path("analysis"))
    p.add_argument("--tree", type=Path, default=Path("Drosophila_tree.nwk"))
    p.add_argument("--focal", default="Dmel",
                   help="focal species leaf, drawn as its own reference row")
    p.add_argument("--order",
                   choices=["elongation", "complete", "cluster", "name",
                            "nested"],
                   default="elongation",
                   help="how to order the gene columns. 'nested' puts the "
                        "monophyletic profiles first, shallowest event first, "
                        "which is the order select_nested.py ranks on")
    p.add_argument("--width", type=float, default=24.0)
    p.add_argument("--suffix", default="",
                   help="appended to the output filenames, to keep two runs "
                        "side by side (e.g. --suffix _flankhigh)")
    return p.parse_args()


# ---------------------------------------------------------------- tree

def tree_layout(newick_path: Path, cladogram: bool = True):
    """Return (leaf_order, segments, xmax) for drawing the phylogeny.

    leaf_order  leaf names, bottom row first (matches origin="lower")
    segments    list of (x0, y0, x1, y1) line segments in tree coordinates
    """
    tree = Phylo.read(str(newick_path), "newick")
    tree.ladderize()

    depths = tree.depths(unit_branch_lengths=cladogram)
    if not cladogram and set(depths.values()) == {0}:
        print("note: the tree has no branch lengths, drawing a cladogram",
              file=sys.stderr)
        depths = tree.depths(unit_branch_lengths=True)
    # pull every tip out to the same depth so the leaves line up against the
    # matrix rather than stopping at ragged distances from it
    dmax = max(depths.values())
    for clade in tree.get_terminals():
        depths[clade] = dmax

    ypos = {}
    for i, clade in enumerate(reversed(tree.get_terminals())):
        ypos[id(clade)] = float(i)

    def assign(clade):
        if clade.is_terminal():
            return ypos[id(clade)]
        ys = [assign(c) for c in clade.clades]
        y = sum(ys) / len(ys)
        ypos[id(clade)] = y
        return y

    assign(tree.root)

    segments = []

    def walk(clade):
        x = depths[clade]
        y = ypos[id(clade)]
        for child in clade.clades:
            cx, cy = depths[child], ypos[id(child)]
            segments.append((x, y, x, cy))    # vertical connector
            segments.append((x, cy, cx, cy))  # horizontal branch
            walk(child)

    walk(tree.root)
    leaf_order = [c.name for c in reversed(tree.get_terminals())]
    return leaf_order, segments, max(depths.values())


# ---------------------------------------------------------------- data

def load_states(path: Path, side: str, pairs: Path | None):
    """Read the canonical states table produced by build_candidates.py.

    All the genomic work — category, verdict, conservative aggregation per
    species — was already done there.
    """
    st = pl.read_csv(path, separator="\t", infer_schema_length=10000)
    for need in ("qseqid", "species", "state"):
        if need not in st.columns:
            sys.exit(f"{path} has no '{need}' column; is it from "
                     f"build_candidates.py?")
    if "side" in st.columns:
        st = st.filter(pl.col("side") == side)
    if st.is_empty():
        sys.exit(f"no rows for side '{side}' in {path}")

    # label the columns by locus rather than by transcript when possible
    if pairs is not None and pairs.exists():
        p = pl.read_csv(pairs, separator="\t", infer_schema_length=10000)
        cols = [c for c in ("qseqid", "gene_name", "gene_id") if c in p.columns]
        if len(cols) >= 2:
            label = (p.select(cols)
                      .unique(subset=["qseqid"], keep="first")
                      .with_columns(
                          pl.coalesce([pl.col(c) for c in cols[1:]] +
                                      [pl.col("qseqid")]).alias("label")))
            st = st.join(label.select(["qseqid", "label"]), on="qseqid",
                         how="left")
    if "label" not in st.columns:
        st = st.with_columns(pl.col("qseqid").alias("label"))
    st = st.with_columns(pl.col("label").fill_null(pl.col("qseqid")))

    return sorted(st["qseqid"].unique().to_list()), st


def state_matrix(genes, species_order, st, states, focal):
    """Dense (species x genes) integer matrix of state codes.

    The focal row is filled with its own code: it has no row in states.tsv
    because build_candidates.py describes the neighbours, and painting it
    'absent' would read as 'no homolog' for the very species the extension
    comes from.
    """
    code_of = {name: v[0] for name, v in states.items()}
    unknown = sorted(set(st["state"].drop_nulls().to_list()) - set(code_of))
    if unknown:
        print(f"WARNING: states absent from this side's palette, drawn as "
              f"'absent': {unknown}", file=sys.stderr)

    grid = (
        pl.DataFrame({"qseqid": [g for g in genes for _ in species_order],
                      "species": species_order * len(genes)})
        .join(st.select(["qseqid", "species", "state"]),
              on=["qseqid", "species"], how="left")
        .with_columns(
            pl.col("state").fill_null("absent")
              .replace_strict(code_of, default=code_of["absent"]).alias("code"))
    )

    wide = grid.pivot(on="qseqid", index="species", values="code")
    rank = {s: i for i, s in enumerate(species_order)}
    wide = (wide.with_columns(pl.col("species").replace_strict(rank).alias("_r"))
                .sort("_r").drop("_r"))
    mat = wide.select(genes).to_numpy()

    if focal in rank:
        mat[rank[focal], :] = code_of["focal"]
    else:
        print(f"note: focal species '{focal}' is not a leaf of the tree",
              file=sys.stderr)
    return mat, grid


def order_genes(mat, genes, how, states, nested_path):
    """Column order. Returns an index array over `genes`."""
    if how == "name":
        return np.argsort(genes)

    if how == "cluster":
        try:
            from scipy.cluster.hierarchy import linkage, leaves_list
            from scipy.spatial.distance import pdist
            return leaves_list(linkage(pdist(mat.T, metric="hamming"),
                                       method="average"))
        except ImportError:
            print("scipy unavailable, falling back to --order elongation",
                  file=sys.stderr)
            return order_genes(mat, genes, "elongation", states, nested_path)

    if how == "nested":
        if nested_path is None or not nested_path.exists():
            print(f"missing {nested_path}; run select_nested.py first. "
                  f"Falling back to --order elongation", file=sys.stderr)
            return order_genes(mat, genes, "elongation", states, nested_path)
        nd = pl.read_csv(nested_path, separator="\t", infer_schema_length=10000)
        mono = dict(zip(nd["qseqid"].to_list(), nd["monophyletic"].to_list()))
        depth = dict(zip(nd["qseqid"].to_list(), nd["event_depth"].to_list()))
        n_el = (mat == states["elongation"][0]).sum(axis=0)
        # monophyletic first, then the shallowest event, then breadth of trace.
        # Genes absent from the nested table were filtered out there; they go
        # last rather than being dropped, since this script does not select.
        key = np.array([
            (0 if mono.get(g) else 1,
             depth.get(g) if depth.get(g) is not None else 99,
             -int(n_el[i]))
            for i, g in enumerate(genes)],
            dtype=[("m", int), ("d", int), ("e", int)])
        return np.argsort(key, order=("m", "d", "e"))

    if how == "complete":
        return np.argsort(-(mat == states["complete"][0]).sum(axis=0))

    # default: breadth of recovered trace first, ties broken by fewest
    # 'complete' species, which are the ones arguing against the candidate
    n_el = (mat == states["elongation"][0]).sum(axis=0)
    n_cp = (mat == states["complete"][0]).sum(axis=0)
    return np.lexsort((n_cp, -n_el))


# ---------------------------------------------------------------- plot

def draw(mat, genes, species_order, segments, xmax, focal, outpath, width,
         states, side):
    n_sp, n_g = mat.shape

    ordered = sorted(states.values(), key=lambda v: v[0])
    codes = [c for c, _, _ in ordered]
    cmap = ListedColormap([c for _, c, _ in ordered])
    norm = BoundaryNorm([c - 0.5 for c in codes] + [max(codes) + 0.5], cmap.N)

    height = max(4.5, 0.30 * n_sp + 2.2)
    fig, (axt, axm) = plt.subplots(
        1, 2, figsize=(width, height),
        gridspec_kw={"width_ratios": [1, 9], "wspace": 0.22})

    # --- tree: root on the left, tips on the right pointing at the matrix
    for x0, y0, x1, y1 in segments:
        axt.plot([x0, x1], [y0, y1], color="#333333", lw=0.9,
                 solid_capstyle="round")
    axt.set_xlim(-0.02 * xmax, xmax * 1.02)
    axt.set_ylim(-0.5, len(species_order) - 0.5)
    axt.axis("off")

    # --- matrix
    axm.imshow(mat, aspect="auto", cmap=cmap, norm=norm,
               interpolation="nearest", origin="lower",
               extent=(-0.5, n_g - 0.5, -0.5, n_sp - 0.5))
    axm.set_yticks(range(n_sp))
    axm.set_yticklabels(
        [s + ("  (focale)" if s == focal else "") for s in species_order],
        fontsize=7)
    axm.tick_params(axis="y", length=0)
    axm.set_xticks([])
    axm.set_xlabel(f"{n_g} gènes candidats", fontsize=9)
    for y in range(n_sp + 1):
        axm.axhline(y - 0.5, color="white", lw=0.4)

    handles = [Patch(facecolor=c, edgecolor="#999999", label=lab)
               for _, c, lab in ordered]
    axm.legend(handles=handles, loc="upper center",
               bbox_to_anchor=(0.5, -0.04), ncol=min(len(handles), 4),
               fontsize=8, frameon=False)

    fig.suptitle(f"Candidats {side.upper()} : état de chaque espèce voisine, "
                 "ordonnées selon l'arbre des espèces", fontsize=11, y=0.98)
    fig.savefig(outpath, bbox_inches="tight", dpi=200)
    plt.close(fig)


# ---------------------------------------------------------------- main

def main():
    args = parse_args()
    states = STATES_CTER if args.side == "cter" else STATES_NTER
    nested_path = args.nested or (args.states.parent /
                                  f"{args.side}_nested.tsv")

    if not args.states.exists():
        sys.exit(f"missing {args.states}; run build_candidates.py first")
    if not args.tree.exists():
        sys.exit(f"missing tree: {args.tree}")

    args.outdir.mkdir(parents=True, exist_ok=True)

    leaf_order, segments, xmax = tree_layout(args.tree, args.cladogram)
    genes, st = load_states(args.states, args.side, args.pairs)

    if args.genes_from:
        if not args.genes_from.exists():
            sys.exit(f"missing {args.genes_from}")
        wanted = [l.strip() for l in open(args.genes_from) if l.strip()]
        keep = [g for g in genes if g in set(wanted)]
        unseen = sorted(set(wanted) - set(genes))
        if unseen:
            print(f"note: {len(unseen)} id(s) from {args.genes_from} absent "
                  f"from states.tsv", file=sys.stderr)
        if not keep:
            sys.exit(f"none of the ids in {args.genes_from} is in states.tsv")
        genes = keep
        st = st.filter(pl.col("qseqid").is_in(genes))

    missing = set(st["species"].drop_nulls().to_list()) - set(leaf_order)
    if missing:
        print(f"WARNING: species absent from the tree, not drawn: "
              f"{sorted(missing)}", file=sys.stderr)

    mat, grid = state_matrix(genes, leaf_order, st, states, args.focal)

    idx = order_genes(mat, genes, args.order, states, nested_path)
    mat = mat[:, idx]
    genes_ord = [genes[i] for i in idx]

    out_pdf = args.outdir / f"{args.side}_tree_matrix{args.suffix}.pdf"
    draw(mat, genes_ord, leaf_order, segments, xmax, args.focal,
         out_pdf, args.width, states, args.side)

    # ---- what the figure contains, in numbers
    rev = {v[0]: k for k, v in states.items()}
    counts = (grid.group_by("code").len().sort("code")
              .with_columns(pl.col("code").replace_strict(rev).alias("state")))
    total = grid.height
    print(f"\n{args.side.upper()}  {len(genes_ord)} genes x "
          f"{len(leaf_order) - 1} neighbour species\n")
    print(counts.select(["state", "len"]).with_columns(
        (100 * pl.col("len") / total).round(1).alias("pct")))

    # the cost of treating coding_flank as decisive, per gene
    cf = (grid.filter(pl.col("state") == "coding_flank")
              .group_by("qseqid").len().rename({"len": "n_coding_flank"}))
    if not cf.is_empty():
        print(f"\ncoding_flank: {cf.height}/{len(genes_ord)} genes affected, "
              f"median {int(cf['n_coding_flank'].median())} species per "
              f"affected gene, max {int(cf['n_coding_flank'].max())}")

    # gene order as drawn, so a column can be traced back to its locus
    labels = dict(zip(st["qseqid"].to_list(), st["label"].to_list()))
    cf_map = dict(zip(cf["qseqid"].to_list(),
                      cf["n_coding_flank"].to_list())) if not cf.is_empty() else {}
    out_tsv = args.outdir / f"{args.side}_column_order{args.suffix}.tsv"
    pl.DataFrame({"column": list(range(1, len(genes_ord) + 1)),
                  "qseqid": genes_ord,
                  "label": [labels.get(g, g) for g in genes_ord],
                  "n_coding_flank": [cf_map.get(g, 0) for g in genes_ord]}
                 ).write_csv(out_tsv, separator="\t")

    print(f"\nfigure : {out_pdf}")
    print(f"columns: {out_tsv}")


if __name__ == "__main__":
    main()