#!/usr/bin/env python3
"""
parse_big_table.py - Parse ssearch results for big send C-terminal candidates.

This script reads nucleotide alignment results and calculates the recovered
elongation for each query-subject pair.

Usage:
    python3 parse_big_table.py <ssearch_big_table.tsv> <qends.tsv>

Output:
    {qseqid}_final_res.tsv
"""

import sys
import os
import polars as pl


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <ssearch_big_table> <qends>", file=sys.stderr)
        sys.exit(1)

    ssearch_big_table = sys.argv[1]
    qends_path = sys.argv[2]

    schema = {
            "qseqid": pl.Utf8,
            "sseqid": pl.Utf8,
            "qend_mapped": pl.Int32,
            "qend_nuc": pl.Int32,
            "qend_difference": pl.Int32,
            "qend": pl.Int32,
            "recovered_elongation": pl.Float32,
            "raw_recovered_elongation": pl.Int32,
            "category": pl.Utf8,
            "qend_prot": pl.Int32
        }

    
    columns = [
        "qseqid",
        "sseqid",
        "qend_mapped",
        "qend_nuc",
        "qend_difference",
        "qend",
        "recovered_elongation",
        "raw_recovered_elongation",
        "category",
        "qend_prot"
    ]


    qseqid = os.path.basename(ssearch_big_table).split("_big")[0]

    if os.path.getsize(ssearch_big_table) == 0 or sum(1 for _ in open(ssearch_big_table)) <= 1:
        pl.DataFrame(schema=schema).write_csv(
            f"{qseqid}_final_res.tsv", separator="\t", include_header=True)
        sys.exit(0)

    table = (
        pl.read_csv(
            ssearch_big_table,
            separator="\t",
            has_header=True,
            schema_overrides={"evalue": pl.Float32}
        )
        .sort(["qseqid", "sseqid", "bitscore"], descending=[False, False, True])
        .unique(subset=["qseqid", "sseqid"], keep="first", maintain_order=True)
        .select(["qseqid", "sseqid", "qend_nuc"])
    )

    if table.height == 0:
        pl.DataFrame(schema=schema).write_csv(
            f"{qseqid}_final_res.tsv", separator="\t", include_header=True)
        sys.exit(0)

    qends_mapped = pl.read_csv(qends_path, separator="\t", has_header=True)

    # ssearch only reports subjects that produced a nucleotide alignment, so its
    # subjects are a subset of qends (which lists every subject). Only the reverse
    # — ssearch carrying a subject absent from qends — signals a real plumbing bug.
    if table.select(pl.col("sseqid")).n_unique() > qends_mapped.select(pl.col("sseqid")).n_unique():
        raise ValueError(
            "The ssearch output has more unique subjects than the qends_mapped file "
            "(ssearch subjects should be a subset of qends subjects)"
        )

    

    (
        table
        .join(qends_mapped, on="sseqid", how="inner")
        .with_columns(
            qend_difference=pl.col("qend_mapped") - pl.col("qend_nuc"),
            qend_prot=pl.col("qend"),
            recovered_elongation=(
                (pl.col("qend_nuc") - ((pl.col("qend") - 1) * 3 + 1))
                / (pl.col("qlend") * 3)
            ),
            raw_recovered_elongation=(pl.col("qend_nuc") - ((pl.col("qend") - 1) * 3 + 1)),
            category=pl.lit("big")
        )
        .select(columns)
        .cast(schema)
        .write_csv(f"{qseqid}_final_res.tsv", separator="\t", include_header=True)
    )


if __name__ == "__main__":
    main()
