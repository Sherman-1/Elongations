#!/usr/bin/env python3
"""
plot_all_matrices.py - One figure gathering the gene x species matrices of every
dataset, N-terminal and C-terminal side by side.

Each dataset is one row. The species tree is drawn once on the left, the species
names once, and the two matrices share the same row order, so a species reads
straight across both extremities. One legend for the whole figure.

This script draws, it does not compute. The states come from each dataset's own
analysis/tables/states.tsv, and the palettes, layout of the tree and ordering of
the genes are imported from tree_matrix.py, which must sit next to this file.
Whatever you changed there (French labels, colour of the focal row) is picked up
here.

Datasets are given as NAME=DIR:TREE:FOCAL, the tree path being relative to DIR
unless absolute. NAME is what is printed above the row.

Usage:
    python3 plot_all_matrices.py \\
        --dataset "D. melanogaster=.:Drosophila_tree.nwk:Dmel" \\
        --dataset "M. musculus=../Elongations_Mmus:Mammal_tree.nwk:Mmus" \\
        --dataset "S. cerevisiae=../Elongations_Scer:Saccharomyces_tree.nwk:Scer" \\
        --abbrev --out analysis/all_matrices.pdf
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.patches import Patch

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
try:
    import tree_matrix as tm
except ImportError:
    sys.exit("tree_matrix.py must sit next to this script")


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dataset", action="append", required=True,
                   metavar="NAME=DIR:TREE:FOCAL")
    p.add_argument("--order", choices=["elongation", "complete", "cluster",
                                       "name"],
                   default="elongation",
                   help="gene column order, as in tree_matrix.py")
    p.add_argument("--abbrev", action="store_true",
                   help="write Drosophila_biarmipes as D. biarmipes. Only long "
                        "genus prefixes are shortened, so Spar_NCBI stays as is")
    p.add_argument("--width", type=float, default=7.0,
                   help="figure width in inches; 7 fits an A4 text block")
    p.add_argument("--row-height", type=float, default=0.13,
                   help="inches per species row")
    p.add_argument("--fontsize", type=float, default=5.5)
    p.add_argument("--out", type=Path, default=Path("analysis/all_matrices.pdf"))
    return p.parse_args()


def parse_dataset(spec):
    if "=" not in spec:
        sys.exit(f"--dataset expects NAME=DIR:TREE:FOCAL, got '{spec}'")
    name, rest = spec.split("=", 1)
    parts = rest.split(":")
    if len(parts) != 3:
        sys.exit(f"--dataset expects NAME=DIR:TREE:FOCAL, got '{spec}'")
    d, tree, focal = Path(parts[0]), Path(parts[1]), parts[2]
    if not tree.is_absolute():
        tree = d / tree
    for f in (d / "analysis" / "tables" / "states.tsv", tree):
        if not f.exists():
            sys.exit(f"missing {f}")
    return name, d, tree, focal


def abbreviate(s):
    if "_" in s:
        genus, rest = s.split("_", 1)
        if len(genus) > 4 and genus[0].isupper():
            return f"{genus[0]}. {rest.replace('_', ' ')}"
    return s


def panel(ddir, tree, focal, side, order):
    states = tm.STATES_CTER if side == "cter" else tm.STATES_NTER
    leaf_order, segments, xmax = tm.tree_layout(tree, True)
    genes, st = tm.load_states(ddir / "analysis" / "tables" / "states.tsv",
                               side, None)
    mat, _ = tm.state_matrix(genes, leaf_order, st, states, focal)
    idx = tm.order_genes(mat, genes, order, states, None)
    return {"mat": mat[:, idx], "leaves": leaf_order, "segments": segments,
            "xmax": xmax, "states": states}


def draw_matrix(ax, mat, states):
    ordered = sorted(states.values(), key=lambda v: v[0])
    codes = [c for c, _, _ in ordered]
    cmap = ListedColormap([c for _, c, _ in ordered])
    norm = BoundaryNorm([c - 0.5 for c in codes] + [max(codes) + 0.5], cmap.N)
    n_sp, n_g = mat.shape
    # 'none' embeds the matrix at full resolution in a PDF: with a thousand
    # columns squeezed into a few centimetres, any resampling here would drop
    # the rare red columns, which are precisely the ones worth seeing
    ax.imshow(mat, aspect="auto", cmap=cmap, norm=norm, interpolation="none",
              origin="lower", extent=(-0.5, n_g - 0.5, -0.5, n_sp - 0.5))
    for y in range(n_sp + 1):
        ax.axhline(y - 0.5, color="white", lw=0.3)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_linewidth(0.4)


def main():
    args = parse_args()
    fs = args.fontsize
    datasets = [parse_dataset(s) for s in args.dataset]

    rows = []
    for name, d, tree, focal in datasets:
        print(f"{name}: {d}", file=sys.stderr)
        nter = panel(d, tree, focal, "nter", args.order)
        cter = panel(d, tree, focal, "cter", args.order)
        rows.append((name, focal, nter, cter))

    heights = [len(n["leaves"]) for _, _, n, _ in rows]
    fig_h = args.row_height * sum(heights) + 0.45 * len(rows) + 0.9
    fig = plt.figure(figsize=(args.width, fig_h))
    outer = GridSpec(len(rows), 1, figure=fig, height_ratios=heights,
                     hspace=0.45, left=0.01, right=0.99, top=0.96,
                     bottom=0.75 / fig_h)

    for r, (name, focal, nter, cter) in enumerate(rows):
        inner = GridSpecFromSubplotSpec(
            1, 4, subplot_spec=outer[r],
            width_ratios=[0.9, 1.8, 6, 6], wspace=0.03)
        ax_tree = fig.add_subplot(inner[0])
        ax_lab = fig.add_subplot(inner[1])
        ax_n = fig.add_subplot(inner[2])
        ax_c = fig.add_subplot(inner[3])

        leaves = nter["leaves"]
        n = len(leaves)

        for x0, y0, x1, y1 in nter["segments"]:
            ax_tree.plot([x0, x1], [y0, y1], color="#333333", lw=0.5,
                         solid_capstyle="round")
        ax_tree.set_xlim(-0.02 * nter["xmax"], nter["xmax"] * 1.02)
        ax_tree.set_ylim(-0.5, n - 0.5)
        ax_tree.axis("off")
        ax_tree.set_title(name, loc="left", fontsize=fs + 2, style="italic",
                          pad=3)

        ax_lab.set_xlim(0, 1)
        ax_lab.set_ylim(-0.5, n - 0.5)
        ax_lab.axis("off")
        for y, sp in enumerate(leaves):
            label = abbreviate(sp) if args.abbrev else sp
            if sp == focal:
                label += " (focale)"
            ax_lab.text(0.03, y, label, va="center", ha="left", fontsize=fs,
                        fontweight="bold" if sp == focal else "normal")

        draw_matrix(ax_n, nter["mat"], nter["states"])
        draw_matrix(ax_c, cter["mat"], cter["states"])
        ax_n.set_xlabel(f"{nter['mat'].shape[1]} gènes candidats",
                        fontsize=fs, labelpad=2)
        ax_c.set_xlabel(f"{cter['mat'].shape[1]} gènes candidats",
                        fontsize=fs, labelpad=2)
        if r == 0:
            ax_n.set_title("Extensions N-terminales", fontsize=fs + 2, pad=8)
            ax_c.set_title("Extensions C-terminales", fontsize=fs + 2, pad=8)

    # one legend: the N-ter palette is a superset of the C-ter one, with the
    # same colour for every shared state
    seen, handles = set(), []
    for pal in (tm.STATES_NTER, tm.STATES_CTER):
        for key, (_, colour, label) in sorted(pal.items(),
                                              key=lambda kv: kv[1][0]):
            if key in seen:
                continue
            seen.add(key)
            handles.append(Patch(facecolor=colour, edgecolor="#999999",
                                 linewidth=0.4, label=label))
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=fs + 0.5,
               frameon=False, bbox_to_anchor=(0.5, 0.0),
               handlelength=1.4, columnspacing=1.5)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, bbox_inches="tight", dpi=600)
    plt.close(fig)
    print(f"figure: {args.out}  ({args.width:.1f} x {fig_h:.1f} in)")


if __name__ == "__main__":
    main()
