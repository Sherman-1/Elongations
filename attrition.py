#!/usr/bin/env python3
"""
attrition.py - How many candidates survive each step, and what each filter costs.

Every criterion in the pipeline removes something. Stating how much is what
justifies its existence: a filter that removes nothing is decoration, one that
removes almost everything needs an argument. The table produced here is also the
only way to compare two datasets, since the raw counts scale with the size of the
focal proteome.

Steps, and where each count comes from:

  1 protéines avec homologue eucaryote   rows of test/statistics
  2 candidats après le crible à 5%       test/{side}_candidates
  3 candidats avec homologue local       distinct qseqid in output/{side}/full.blast
  4 candidats avec une espèce retenue    distinct qseqid in smallFinal or bigFinal
  5 candidats en catégorie small         distinct qseqid in smallFinal
  6 candidats avec une trace             qseqid with a verdict of elongation
  7 loci après déduplication             genes.tsv, key (gene_id, ext_aa)
  8 loci soutenus par N espèces          see below
  9 et extension dans l'exon terminal    ext_in_terminal_exon

Steps 1 to 6 count TRANSCRIPTS, steps 7 to 9 count LOCI: the drop between 6 and
7 is isoform redundancy, not selection. It is reported as such rather than
folded into the funnel.

Step 8 can be counted two ways, and they are NOT interchangeable:

  n_species_elongation   species whose aggregated verdict is 'elongation', i.e.
                         geometry AND nothing more conservative overriding it.
                         A species carrying an annotation_doubt does not count
  n_species_contiguous   species with a contiguous HSP over the extension,
                         geometry alone, regardless of the verdict

genes.tsv is sorted on the second while the funnel has historically filtered on
the first, so a locus can sit at the top of one and below the bar of the other.
Both are printed; --support-column decides which one drives steps 8 and 9.

Usage:
    python3 attrition.py --side nter
    python3 attrition.py --side nter --dataset Dmel=. --dataset Mmus=../Elongations_Mmus
"""

import argparse
import sys
from pathlib import Path

import polars as pl
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--side", choices=["nter", "cter"], default="nter")
    p.add_argument("--dataset", action="append", default=None,
                   metavar="NAME=PATH",
                   help="repeatable; defaults to a single dataset in the "
                        "current directory")
    p.add_argument("--min-species", type=int, default=5,
                   help="species supporting a locus at the last step")
    p.add_argument("--support-column",
                   choices=["n_species_elongation", "n_species_contiguous"],
                   default="n_species_elongation",
                   help="which support counter drives steps 8 and 9")
    p.add_argument("--outdir", type=Path, default=Path("analysis"))
    return p.parse_args()


def n_lines(path: Path, header=True):
    if not path.exists():
        return None
    with open(path) as fh:
        n = sum(1 for _ in fh)
    return max(0, n - 1) if header else n


def n_unique(path: Path, col="qseqid", sep="\t"):
    if not path.exists():
        return None
    try:
        return pl.read_csv(path, separator=sep,
                           infer_schema_length=10000)[col].n_unique()
    except Exception as e:                                     # noqa: BLE001
        print(f"  note: could not read {path}: {e}", file=sys.stderr)
        return None


def counts_for(root: Path, side: str, min_species: int, support_col: str):
    """Ordered list of (label, count, unit). None means the file was not found."""
    out = []
    tables = root / "analysis" / "tables"

    out.append(("protéines avec homologue eucaryote",
                n_lines(root / "test" / "statistics"), "transcrits"))
    out.append(("candidats après le crible à 5%",
                n_lines(root / "test" / f"{side}_candidates", header=False),
                "transcrits"))
    out.append(("candidats avec homologue local",
                n_unique(root / "output" / side / "full.blast"), "transcrits"))

    small = root / "output" / side / "smallFinal"
    big = root / "output" / side / "bigFinal"
    kept = None
    a = n_unique(small)
    b = n_unique(big)
    if a is not None or b is not None:
        s = set()
        for f in (small, big):
            if f.exists():
                s |= set(pl.read_csv(f, separator="\t",
                                     infer_schema_length=10000)["qseqid"].to_list())
        kept = len(s)
    out.append(("candidats avec une espèce retenue", kept, "transcrits"))
    out.append(("candidats en catégorie small", a, "transcrits"))

    pairs = tables / "pairs.tsv"
    n_trace = None
    if pairs.exists():
        p = pl.read_csv(pairs, separator="\t", infer_schema_length=10000)
        p = p.filter(pl.col("side") == side)
        if "verdict" in p.columns:
            n_trace = p.filter(pl.col("verdict") == "elongation")["qseqid"].n_unique()
    out.append(("candidats avec une trace", n_trace, "transcrits"))

    genes = tables / "genes.tsv"
    n_loci = n_sup = n_exon = None
    alt = {}
    if genes.exists():
        g = pl.read_csv(genes, separator="\t", infer_schema_length=10000)
        g = g.filter(pl.col("side") == side)
        n_loci = g["gene_id"].drop_nulls().n_unique()
        for col in ("n_species_elongation", "n_species_contiguous"):
            if col in g.columns:
                alt[col] = g.filter(pl.col(col) >= min_species).height
        if support_col in g.columns:
            sup = g.filter(pl.col(support_col) >= min_species)
            n_sup = sup.height
            if "ext_in_terminal_exon" in g.columns:
                n_exon = sup.filter(
                    pl.col("ext_in_terminal_exon")
                      .cast(pl.Boolean, strict=False).fill_null(False)).height
        else:
            print(f"  note: {support_col} absent from {genes}", file=sys.stderr)
    out.append(("loci après déduplication", n_loci, "loci"))
    out.append((f"loci soutenus par >= {min_species} espèces", n_sup, "loci"))
    out.append(("et extension dans l'exon terminal", n_exon, "loci"))
    return out, alt


def main():
    args = parse_args()
    datasets = []
    if args.dataset:
        for d in args.dataset:
            if "=" not in d:
                sys.exit(f"--dataset expects NAME=PATH, got '{d}'")
            name, path = d.split("=", 1)
            datasets.append((name, Path(path)))
    else:
        datasets.append((Path.cwd().name, Path(".")))

    tables, alts = {}, {}
    for name, root in datasets:
        if not root.is_dir():
            sys.exit(f"{root} is not a directory")
        tables[name], alts[name] = counts_for(root, args.side, args.min_species,
                                              args.support_column)

    labels = [lab for lab, _, _ in tables[datasets[0][0]]]
    units = [u for _, _, u in tables[datasets[0][0]]]

    # ---- text table
    w = max(len(l) for l in labels) + 2
    head = f"{'étape':<{w}}" + "".join(f"{n:>22}" for n, _ in datasets)
    print(f"\n{args.side.upper()}   (support: {args.support_column})\n")
    print(head)
    print("-" * len(head))

    rows_out = {"etape": labels, "unite": units}
    for name, _ in datasets:
        rows_out[name] = [c for _, c, _ in tables[name]]

    prev = {name: None for name, _ in datasets}
    for i, lab in enumerate(labels):
        line = f"{lab:<{w}}"
        for name, _ in datasets:
            c = tables[name][i][1]
            if c is None:
                line += f"{'-':>22}"
            else:
                # percentage of the previous step, within the same unit
                same_unit = i > 0 and units[i] == units[i - 1]
                p = prev[name]
                pct = (f"  ({100*c/p:.0f}%)"
                       if same_unit and p not in (None, 0) else "")
                line += f"{str(c) + pct:>22}"
                prev[name] = c
        print(line)

    # ---- the two support counters side by side, because they disagree
    if any(alts.values()):
        print(f"\nloci >= {args.min_species} espèces, selon le compteur:")
        for col in ("n_species_elongation", "n_species_contiguous"):
            vals = "  ".join(f"{name}={alts[name].get(col, '-')}"
                             for name, _ in datasets)
            mark = " <- utilisé" if col == args.support_column else ""
            print(f"  {col:<24} {vals}{mark}")

    args.outdir.mkdir(parents=True, exist_ok=True)
    out_tsv = args.outdir / f"{args.side}_attrition.tsv"
    pl.DataFrame(rows_out).write_csv(out_tsv, separator="\t")

    # ---- figure: one bar per step, log scale because the range spans 4 orders
    fig, ax = plt.subplots(figsize=(9, 0.5 * len(labels) + 2))
    n_ds = len(datasets)
    h = 0.8 / n_ds
    colors = ["#2f4b7c", "#c3423f", "#3f7d5c", "#e8a33d"]
    for k, (name, _) in enumerate(datasets):
        vals = [tables[name][i][1] for i in range(len(labels))]
        ypos = [i + (k - (n_ds - 1) / 2) * h for i in range(len(labels))]
        vv = [v if v else 0.5 for v in vals]      # log scale needs > 0
        ax.barh(ypos, vv, height=h * 0.9, label=name,
                color=colors[k % len(colors)])
        for y, v in zip(ypos, vals):
            if v is not None:
                ax.text(max(v, 0.5) * 1.15, y, str(v), va="center", fontsize=7)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlabel("nombre (échelle log)", fontsize=9)
    if n_ds > 1:
        ax.legend(fontsize=8, frameon=False)
    ax.set_title(f"{args.side.upper()} — attrition des candidats", fontsize=10)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    out_pdf = args.outdir / f"{args.side}_attrition.pdf"
    fig.savefig(out_pdf, bbox_inches="tight", dpi=200)
    plt.close(fig)

    print(f"\ntableau: {out_tsv}")
    print(f"figure : {out_pdf}")
    missing = [labels[i] for i in range(len(labels))
               if all(tables[n][i][1] is None for n, _ in datasets)]
    if missing:
        print(f"\nnon trouvé: {', '.join(missing)}", file=sys.stderr)


if __name__ == "__main__":
    main()