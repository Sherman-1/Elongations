#!/usr/bin/env python3
"""
screen_stats.py - What the eukaryote-scale screen sees, before any neighbour is
looked at. Produces the numbers and figures for the first results section.

Reads test/statistics, one row per focal protein with a eukaryotic homolog:

  nb_unique_species    species with at least one homolog passing the homology
                       filters. This is the DENOMINATOR of the screen
  nter_q_15            among them, species with an alignment reaching within 15
                       residues of the N-terminus; cter_q_15 likewise
  percent_*_q_15       the same as a percentage. A protein is a candidate when
                       this stays at or below the threshold
  qstart_median        median first aligned residue across homologs — the
                       N-terminal extension length as the screen sees it
  cter_median          the same at the other end

Three things are reported.

DENOMINATOR. A percentage threshold does not mean the same thing at both ends
of the species-count distribution. At 16 species, 5% allows no covering species
at all; at 2000, it allows a hundred. The screen is therefore stricter on poorly
represented proteins than on well represented ones, and the distribution of
nb_unique_species is what makes that visible.

OVERLAP. Proteins candidate at BOTH ends are counted, and compared to what
independence would give. An excess says the screen partly captures a property of
the protein rather than an event at one extremity.

LENGTHS. Extension lengths at this scale are relational to NR, not to the
neighbour set analysed later; the two are not the same quantity and the figure
labels say so.

Usage:
    python3 screen_stats.py
    python3 screen_stats.py --dataset Dmel=. --dataset Mmus=../Elongations_Mmus
"""

import argparse
import sys
from pathlib import Path

import polars as pl
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

COLORS = ["#2f4b7c", "#c3423f", "#3f7d5c", "#e8a33d"]


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dataset", action="append", default=None,
                   metavar="NAME=PATH",
                   help="repeatable; defaults to the current directory")
    p.add_argument("--threshold", type=float, default=5.0,
                   help="percentage of covering species defining a candidate")
    p.add_argument("--outdir", type=Path, default=Path("analysis"))
    return p.parse_args()


def load(root: Path):
    f = root / "test" / "statistics"
    if not f.exists():
        sys.exit(f"missing {f}")
    df = pl.read_csv(f, separator="\t", infer_schema_length=10000)
    needed = ["qseqid", "nb_unique_species",
              "percent_nter_q_15", "percent_cter_q_15"]
    missing = [c for c in needed if c not in df.columns]
    if missing:
        sys.exit(f"{f} lacks {missing}")
    return df


def quartiles(s: pl.Series):
    v = s.drop_nulls()
    return (int(v.min()), int(v.quantile(0.25)), int(v.median()),
            int(v.quantile(0.75)), int(v.max()))


def describe(name, df, thr):
    n = df.height
    nter = df.filter(pl.col("percent_nter_q_15") <= thr)
    cter = df.filter(pl.col("percent_cter_q_15") <= thr)
    sn, sc = set(nter["qseqid"]), set(cter["qseqid"])
    both = sn & sc
    expected = len(sn) * len(sc) / n if n else 0

    q = quartiles(df["nb_unique_species"])
    print(f"\n=== {name} ===")
    print(f"protéines avec homologue eucaryote : {n}")
    print(f"espèces par protéine  min {q[0]}  Q1 {q[1]}  médiane {q[2]}  "
          f"Q3 {q[3]}  max {q[4]}")

    # the threshold is not equally strict along that distribution
    allowed = df.select(
        (pl.col("nb_unique_species") * thr / 100).floor().alias("k"))["k"]
    n_zero = int((allowed < 1).sum())
    print(f"  protéines pour lesquelles {thr:g}% n'autorise AUCUNE espèce "
          f"couvrante : {n_zero} ({100*n_zero/n:.1f}%)")

    print(f"candidats N-ter : {len(sn)} ({100*len(sn)/n:.1f}%)")
    print(f"candidats C-ter : {len(sc)} ({100*len(sc)/n:.1f}%)")
    print(f"  N-ter seul {len(sn - sc)}   C-ter seul {len(sc - sn)}   "
          f"les deux {len(both)}")
    print(f"  attendu sous indépendance : {expected:.0f}  "
          f"→ rapport {len(both)/expected:.1f}x" if expected else "")

    for side, sub, col in (("N-ter", nter, "qstart_median"),
                           ("C-ter", cter, "cter_median")):
        if col in sub.columns and sub.height:
            qq = quartiles(sub[col])
            print(f"longueur d'extension {side} (échelle NR, résidus) : "
                  f"médiane {qq[2]}  Q1 {qq[1]}  Q3 {qq[3]}  max {qq[4]}")

    return {
        "dataset": name, "n_proteins": n,
        "species_min": q[0], "species_q1": q[1], "species_median": q[2],
        "species_q3": q[3], "species_max": q[4],
        "n_zero_allowed": n_zero,
        "n_nter": len(sn), "n_cter": len(sc),
        "n_both": len(both), "n_both_expected": round(expected, 1),
        "pct_nter": round(100 * len(sn) / n, 2),
        "pct_cter": round(100 * len(sc) / n, 2),
    }, nter, cter


def figure(datasets, frames, cands, thr, out):
    n = len(datasets)
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))

    # (a) denominator
    ax = axes[0, 0]
    for k, (name, _) in enumerate(datasets):
        v = frames[name]["nb_unique_species"].drop_nulls().to_numpy()
        v = v[v > 0]
        ax.hist(v, bins=np.logspace(0, np.log10(max(v.max(), 10)), 40),
                histtype="step", lw=1.6, color=COLORS[k % 4], label=name)
    ax.set_xscale("log")
    ax.set_xlabel("espèces eucaryotes par protéine", fontsize=9)
    ax.set_ylabel("protéines", fontsize=9)
    ax.set_title("a. dénominateur du crible", fontsize=10, loc="left")
    ax.legend(fontsize=8, frameon=False)

    # (b) the screened quantity, both ends
    ax = axes[0, 1]
    for k, (name, _) in enumerate(datasets):
        d = frames[name]
        for col, ls in (("percent_nter_q_15", "-"), ("percent_cter_q_15", "--")):
            ax.hist(d[col].drop_nulls().to_numpy(), bins=50, histtype="step",
                    lw=1.4, ls=ls, color=COLORS[k % 4],
                    label=f"{name} {'N' if 'nter' in col else 'C'}-ter")
    ax.axvline(thr, color="#888888", lw=1, ls=":")
    ax.set_xlabel("% d'espèces couvrant l'extrémité", fontsize=9)
    ax.set_ylabel("protéines", fontsize=9)
    ax.set_yscale("log")
    ax.set_title(f"b. critère du crible (seuil {thr:g}%)", fontsize=10,
                 loc="left")
    ax.legend(fontsize=7, frameon=False)

    # (c) threshold strictness against the denominator
    ax = axes[1, 0]
    for k, (name, _) in enumerate(datasets):
        d = frames[name]
        ax.scatter(d["nb_unique_species"].to_numpy(),
                   d["percent_nter_q_15"].to_numpy(),
                   s=2, alpha=0.15, color=COLORS[k % 4], edgecolors="none",
                   label=name)
    ax.axhline(thr, color="#888888", lw=1, ls=":")
    ax.set_xscale("log")
    ax.set_xlabel("espèces eucaryotes par protéine", fontsize=9)
    ax.set_ylabel("% couvrant le N-ter", fontsize=9)
    ax.set_title("c. le seuil n'est pas également strict", fontsize=10,
                 loc="left")
    ax.legend(fontsize=8, frameon=False, markerscale=4)

    # (d) extension lengths of the candidates
    ax = axes[1, 1]
    drawn = False
    for k, (name, _) in enumerate(datasets):
        for side, col, ls in (("N", "qstart_median", "-"),
                              ("C", "cter_median", "--")):
            sub = cands[(name, side)]
            if col not in sub.columns or not sub.height:
                continue
            v = sub[col].drop_nulls().to_numpy()
            v = v[v > 0]
            if not len(v):
                continue
            ax.hist(v, bins=np.logspace(0, np.log10(max(v.max(), 10)), 40),
                    histtype="step", lw=1.4, ls=ls, color=COLORS[k % 4],
                    label=f"{name} {side}-ter")
            drawn = True
    if drawn:
        ax.set_xscale("log")
        ax.legend(fontsize=7, frameon=False)
    ax.set_xlabel("longueur d'extension, échelle NR (résidus)", fontsize=9)
    ax.set_ylabel("candidats", fontsize=9)
    ax.set_title("d. longueurs au stade du crible", fontsize=10, loc="left")

    for a in axes.ravel():
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight", dpi=200)
    plt.close(fig)


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

    rows, frames, cands = [], {}, {}
    for name, root in datasets:
        df = load(root)
        frames[name] = df
        r, nter, cter = describe(name, df, args.threshold)
        cands[(name, "N")] = nter
        cands[(name, "C")] = cter
        rows.append(r)

    args.outdir.mkdir(parents=True, exist_ok=True)
    out_tsv = args.outdir / "screen_stats.tsv"
    pl.DataFrame(rows).write_csv(out_tsv, separator="\t")
    out_pdf = args.outdir / "screen_stats.pdf"
    figure(datasets, frames, cands, args.threshold, out_pdf)

    print(f"\ntableau: {out_tsv}")
    print(f"figure : {out_pdf}")
    print("\nles comptes de candidats doivent correspondre à "
          "test/{nter,cter}_candidates ; un écart signale un filtre "
          "supplémentaire appliqué en amont", file=sys.stderr)


if __name__ == "__main__":
    main()