#!/usr/bin/env python3
"""
extract_proteins.py - Focal protein and its neighbours' ANNOTATED proteins, as
a plain multifasta.

This is the simplest comparison there is: the candidate with its extension, and
the proteins the neighbours actually encode. Aligning them shows the extension
facing nothing, which is the observation the whole analysis starts from.

It is NOT the elongated translation. That one shows what the upstream region
WOULD give if it were translated, stops included, and answers a different
question. Use msa_candidates.py for it.

One subject per species by default, chosen on the blastp percent-positives —
the likeliest ortholog, picked independently of whether it supports the
hypothesis.

Usage:
    python3 extract_proteins.py --side nter --genes rna-NM_168110.3
    python3 extract_proteins.py --side nter --from-list top_nter.txt
"""

import argparse
import sys
from pathlib import Path

import polars as pl


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--side", choices=["nter", "cter"], default="nter")
    p.add_argument("--genes", nargs="*", default=None)
    p.add_argument("--from-list", type=Path, default=None,
                   help="file with one gene id per line, e.g. from select_nested.py")
    p.add_argument("--pairs", type=Path,
                   default=Path("analysis/tables/pairs.tsv"))
    p.add_argument("--work", type=Path, default=Path("work"))
    p.add_argument("--proteome", type=Path, default=None,
                   help="neighbour proteome (default: the concatenated.faa "
                        "found under work/)")
    p.add_argument("--outdir", type=Path, default=Path("analysis/proteins"))
    p.add_argument("--all-subjects", action="store_true",
                   help="keep every homolog instead of one per species")
    return p.parse_args()


def read_fasta(path):
    recs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    recs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            elif line:
                buf.append(line)
    if name:
        recs[name] = "".join(buf)
    return recs


def write_fasta(path, records):
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


def find_proteome(work: Path):
    hits = sorted(work.glob("*/*/concatenated.faa"))
    if len(hits) > 1:
        print(f"note: {len(hits)} concatenated.faa under {work}, taking "
              f"{hits[-1]}", file=sys.stderr)
    return hits[-1] if hits else None


def find_query_faa(gene, workdir, work: Path):
    """query.faa lives in the small_elongate_and_align dir.

    build_candidates.py now carries the workdir column through to pairs.tsv, so
    the recorded path is tried first and the glob is only a fallback for tables
    built before that.
    """
    d = work / workdir if workdir else None
    if d and (d / "query.faa").exists():
        return d / "query.faa"
    if d and (d / "elongated_subjects.fna").exists():
        e = (d / "elongated_subjects.fna").resolve().parent
        if (e / "query.faa").exists():
            return e / "query.faa"
    hits = sorted(work.glob(f"*/*/{gene}_subjects_thresholds.tsv"))
    for h in hits:
        if (h.parent / "query.faa").exists():
            return h.parent / "query.faa"
    return None


def main():
    args = parse_args()
    if not args.pairs.exists():
        sys.exit(f"missing {args.pairs}")

    genes = list(args.genes or [])
    if args.from_list:
        genes += [l.strip() for l in open(args.from_list) if l.strip()]
    if not genes:
        sys.exit("give --genes or --from-list")

    prot_path = args.proteome or find_proteome(args.work)
    if prot_path is None or not prot_path.exists():
        sys.exit("neighbour proteome not found; pass --proteome")
    print(f"proteome: {prot_path}", file=sys.stderr)
    proteome = read_fasta(prot_path)

    p = pl.read_csv(args.pairs, separator="\t", infer_schema_length=10000)
    p = p.filter((pl.col("side") == args.side) & (pl.col("category") == "small"))
    if p.is_empty():
        sys.exit(f"no small pair for side '{args.side}' in {args.pairs}")
    has_workdir = "workdir" in p.columns

    args.outdir.mkdir(parents=True, exist_ok=True)
    written = 0

    for gene in genes:
        sub = p.filter(pl.col("qseqid") == gene)
        if sub.is_empty():
            print(f"  {gene}: no small homolog, skipped", file=sys.stderr)
            continue

        if not args.all_subjects and "ppos" in sub.columns:
            # ties broken on the subject id so the pick is reproducible
            sub = (sub.sort(["ppos", "sseqid"], descending=[True, False],
                            nulls_last=True)
                      .unique(subset=["species"], keep="first",
                              maintain_order=True))

        recs = []
        wd = sub["workdir"][0] if has_workdir else None
        qf = find_query_faa(gene, wd, args.work)
        if qf:
            q = read_fasta(qf)
            if q:
                n, s = next(iter(q.items()))
                ext = sub["ext_aa"].max() if "ext_aa" in sub.columns else None
                tag = f"FOCAL_{n}" + (f"_ext{int(ext)}aa" if ext else "")
                recs.append((tag, s))
        else:
            print(f"  {gene}: query.faa not found", file=sys.stderr)

        missing = 0
        for r in sub.iter_rows(named=True):
            sid, sp = r["sseqid"], r.get("species") or "?"
            if sid in proteome:
                recs.append((f"{sp}_{sid}", proteome[sid]))
            else:
                missing += 1

        if len(recs) < 2:
            print(f"  {gene}: fewer than two sequences, skipped", file=sys.stderr)
            continue

        out = args.outdir / f"{gene.replace('/', '_')}_prot.fasta"
        write_fasta(out, recs)
        written += 1
        note = f", {missing} absent(s) du protéome" if missing else ""
        print(f"  {gene}: {len(recs)-1} voisins{note} -> {out}")

    print(f"\n{written} fichier(s) dans {args.outdir}")


if __name__ == "__main__":
    main()