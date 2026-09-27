#!/usr/bin/env python3
"""
select_candidates.py - Pick the candidates whose phylogenetic profile tells the
simplest possible story, straight from states.tsv and the species tree.

WHY NOT THE ORDERED MATRIX. build_ordered_matrix.py flattens the tree into a
line of distance blocks and asks whether the profile is monotonic along it. That
is readable, and it happens to be exact on the Drosophila tree, but it is a
detour: the question is about clades, and clades are read on the topology. This
script skips the matrix entirely. The matrix and tree_matrix.py remain useful
for LOOKING at the data; nothing is selected from them.

It also drops a recoding that was quietly wrong for this purpose. The ordered
matrix sends 'elongation', 'weak_trace' and 'no_trace' to the same value, 0,
because it encodes "is the extension translated". For a reservoir argument that
throws away the observation: a neighbour where the sequence was RECOVERED in the
non-coding flank and one where NOTHING was found are not the same evidence.

THE TEST. Two nested clades.

  clade 1, TRANSLATED   {focal} + the species carrying the extension in an
                        annotated protein ('complete'). If these form a clade,
                        one translation event explains them.

  clade 2, CARRIER      clade 1 + the species where the sequence was recovered
                        in the non-coding flank ('elongation'). If these form a
                        clade, one sequence-origin event explains them.

Clade 1 sits inside clade 2 by construction. A candidate is IDEAL when both are
clades and at least one species is in the reservoir state: the sequence appeared
once, and became translated once, later and closer to the focal species.

The focal species is always in clade 1: it carries the extension by
construction, and has no row in states.tsv.

FOUR ROLES, not seven states. Every state gets exactly one role, listed in
ROLES below and overridable with --neutral:

  coding      'complete'      in both clades
  reservoir   'elongation'    in clade 2 only
  blocking    everything else the pipeline decided ('coding_flank',
                              'weak_trace', 'no_trace', 'annotation_doubt',
                              'not_analysable'). A blocking leaf INSIDE either
                              clade breaks it. This is the strict reading: an
                              undecidable neighbour in the middle of the clade
                              means the story is not simple, whatever the reason
  neutral     no homolog      pruned from the test and counted. An absent
                              homolog is not evidence against anything, and
                              treating it as blocking would fail almost every
                              gene for lack of data rather than for conflict

Monophyly is tested on the informative leaves: a node's descendants minus the
neutral ones must be exactly the set. Pruning is done by intersection rather
than by mutating the tree, which is equivalent and cheaper.

NO LINEARISATION anywhere. Nothing assumes the tree is a comb from the focal
leaf, so this transposes to the mouse tree, where Rnor and Rrat are sisters to
each other at the same distance from Mmus as other species.

Usage:
    python3 select_candidates.py --side nter --focal Dmel --tree Drosophila_tree.nwk
    python3 select_candidates.py --side nter --ideal-only --out-list top_nter.txt
    python3 select_candidates.py --side nter --neutral absent,weak_trace
"""

import argparse
import sys
from pathlib import Path

import polars as pl
from Bio import Phylo

# state -> role. 'absent' is not a state in states.tsv; it is the absence of a
# row, and is handled as neutral below.
ROLES = {
    "complete":         "coding",
    "elongation":       "reservoir",
    "coding_flank":     "blocking",
    "weak_trace":       "blocking",
    "no_trace":         "blocking",
    "annotation_doubt": "blocking",
    "not_analysable":   "blocking",
}
ABSENT = "absent"


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--side", choices=["nter", "cter"], default="nter")
    p.add_argument("--states", type=Path,
                   default=Path("analysis/tables/states.tsv"))
    p.add_argument("--pairs", type=Path,
                   default=Path("analysis/tables/pairs.tsv"),
                   help="only used to label rows by locus")
    p.add_argument("--tree", type=Path, default=Path("Drosophila_tree.nwk"))
    p.add_argument("--focal", default="Dmel")
    p.add_argument("--neutral", default=ABSENT,
                   help=f"comma-separated states to prune from the test instead "
                        f"of letting them break a clade. Default '{ABSENT}' "
                        f"(no homolog). Adding 'weak_trace' is the obvious "
                        f"loosening: it says a trace was seen but under the "
                        f"threshold, which is closer to 'not sure' than to "
                        f"'contradicts'")
    p.add_argument("--min-reservoir", type=int, default=1,
                   help="minimum species in the reservoir state. 1 means "
                        "'at least one neighbour keeps the sequence', which is "
                        "the least that makes the observation exist")
    p.add_argument("--min-informative", type=int, default=1,
                   help="minimum neighbour leaves that are not neutral")
    p.add_argument("--ideal-only", action="store_true",
                   help="keep only the candidates passing both clade tests")
    p.add_argument("--out", type=Path, default=None)
    p.add_argument("--out-list", type=Path, default=None)
    return p.parse_args()


# ---------------------------------------------------------------- tree

def load_tree(path: Path, focal: str):
    tree = Phylo.read(str(path), "newick")
    leaves = {c.name: c for c in tree.get_terminals()}
    if focal not in leaves:
        sys.exit(f"focal species '{focal}' is not a leaf of {path}")
    lineage = [tree.root] + tree.get_path(focal)
    depth_of = {id(c): i for i, c in enumerate(lineage)}
    n = len(lineage) - 1
    return tree, leaves, depth_of, n


def clade_of(tree, leaves, members, neutral, depth_of, n_focal):
    """Is `members` exactly one node's informative descendants?

    Returns (ok, extra, depth). `extra` names the leaves inside the candidate
    clade that are not in `members` and not neutral — the counterexamples.
    `depth` counts edges from the focal leaf up to that node, so 0 means the
    node is the focal leaf itself.
    """
    if not members:
        return False, set(), None
    if len(members) == 1:
        node = leaves[next(iter(members))]
    else:
        node = tree.common_ancestor(*sorted(members))
    covered = {t.name for t in node.get_terminals()} - neutral
    extra = covered - members
    depth = n_focal - depth_of[id(node)] if id(node) in depth_of else None
    return (not extra), extra, depth


# ---------------------------------------------------------------- main

def main():
    args = parse_args()
    if not args.states.exists():
        sys.exit(f"missing {args.states}; run build_candidates.py first")
    if not args.tree.exists():
        sys.exit(f"missing {args.tree}")

    neutral_states = {s.strip() for s in args.neutral.split(",") if s.strip()}
    roles = dict(ROLES)
    for s in neutral_states:
        if s != ABSENT:
            if s not in roles:
                sys.exit(f"--neutral: unknown state '{s}'. Known: "
                         f"{sorted(roles)} and '{ABSENT}'")
            roles[s] = "neutral"

    tree, leaves, depth_of, n_focal = load_tree(args.tree, args.focal)

    st = pl.read_csv(args.states, separator="\t", infer_schema_length=10000)
    for need in ("qseqid", "species", "state"):
        if need not in st.columns:
            sys.exit(f"{args.states} has no '{need}' column")
    if "side" in st.columns:
        st = st.filter(pl.col("side") == args.side)
    if st.is_empty():
        sys.exit(f"no rows for side '{args.side}' in {args.states}")

    unknown = sorted(set(st["state"].drop_nulls().to_list()) - set(roles))
    if unknown:
        sys.exit(f"states with no role defined: {unknown}. Add them to ROLES.")

    off_tree = sorted(set(st["species"].drop_nulls().to_list()) - set(leaves))
    if off_tree:
        print(f"WARNING: species absent from the tree, dropped: {off_tree}",
              file=sys.stderr)
        st = st.filter(~pl.col("species").is_in(off_tree))

    neighbours = set(leaves) - {args.focal}

    labels = {}
    if args.pairs.exists():
        p = pl.read_csv(args.pairs, separator="\t", infer_schema_length=10000)
        if "gene_name" in p.columns:
            u = p.select(["qseqid", "gene_name"]).unique(subset=["qseqid"],
                                                         keep="first")
            labels = dict(zip(u["qseqid"].to_list(), u["gene_name"].to_list()))

    rows = []
    for (gene,), sub in st.group_by(["qseqid"], maintain_order=True):
        by_role = {"coding": set(), "reservoir": set(),
                   "blocking": set(), "neutral": set()}
        seen = set()
        for sp, state in zip(sub["species"].to_list(), sub["state"].to_list()):
            by_role[roles[state]].add(sp)
            seen.add(sp)
        # a neighbour with no row has no homolog: neutral, says nothing
        by_role["neutral"] |= (neighbours - seen)
        if ABSENT in neutral_states:
            neutral = by_role["neutral"]
        else:
            neutral = set()
            by_role["blocking"] |= by_role["neutral"]

        n_info = len(neighbours - neutral)
        if n_info < args.min_informative:
            continue
        if len(by_role["reservoir"]) < args.min_reservoir:
            continue

        coding = by_role["coding"] | {args.focal}
        carrier = coding | by_role["reservoir"]

        ok1, extra1, d1 = clade_of(tree, leaves, coding, neutral,
                                   depth_of, n_focal)
        ok2, extra2, d2 = clade_of(tree, leaves, carrier, neutral,
                                   depth_of, n_focal)

        rows.append({
            "qseqid": gene,
            "label": labels.get(gene, gene),
            "n_coding": len(by_role["coding"]),
            "n_reservoir": len(by_role["reservoir"]),
            "n_blocking": len(by_role["blocking"]),
            "n_neutral": len(neutral),
            "n_informative": n_info,
            "translated_clade": ok1,
            "translated_depth": d1,
            "translated_extra": ";".join(sorted(extra1)),
            "carrier_clade": ok2,
            "carrier_depth": d2,
            "carrier_extra": ";".join(sorted(extra2)),
            "ideal": ok1 and ok2,
        })

    if not rows:
        sys.exit("no gene passes --min-reservoir / --min-informative")

    df = pl.DataFrame(rows)
    if args.ideal_only:
        df = df.filter(pl.col("ideal"))
        if df.is_empty():
            sys.exit("no ideal candidate; loosen --neutral or --min-reservoir")

    # ideal first; then the widest reservoir, which is the observation itself;
    # then the shallowest translation event, a gain on the focal branch being
    # the simplest story of all
    df = df.sort(["ideal", "n_reservoir", "translated_depth"],
                 descending=[True, True, False], nulls_last=True)

    pl.Config.set_tbl_rows(40)
    pl.Config.set_tbl_width_chars(240)
    print(f"\n{args.side.upper()}  focal {args.focal}  "
          f"{len(neighbours)} neighbour leaves")
    print(f"neutral states: {sorted(neutral_states)}")
    print(f"blocking states: "
          f"{sorted(s for s, r in roles.items() if r == 'blocking')}\n")
    print(f"genes evaluated       : {df.height}")
    print(f"translated is a clade : {df.filter(pl.col('translated_clade')).height}")
    print(f"carrier is a clade    : {df.filter(pl.col('carrier_clade')).height}")
    print(f"both (ideal)          : {df.filter(pl.col('ideal')).height}")
    n_term = df.filter(pl.col("ideal") & (pl.col("translated_depth") == 0)).height
    print(f"  translation on the focal branch only: {n_term}\n")
    print(df.head(30))

    out = args.out or args.states.parent / f"{args.side}_candidates.tsv"
    df.write_csv(out, separator="\t")
    print(f"\ntable: {out}")
    if args.out_list:
        args.out_list.write_text("\n".join(df["qseqid"].to_list()) + "\n")
        print(f"list : {args.out_list}")


if __name__ == "__main__":
    main()