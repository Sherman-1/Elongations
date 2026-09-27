#!/usr/bin/env python3
"""
parse_big_table.py - Parse ssearch results for big sstart N-terminal candidates.

This script reads nucleotide alignment results and calculates the recovered
elongation for each query-subject pair.

Usage:
    python3 parse_big_table.py <ssearch_big_table.tsv> <qstarts.tsv>

Output:
    {qseqid}_final_res.tsv
"""

import sys
import os
import polars as pl


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <ssearch_big_table> <qstarts>", file=sys.stderr)
        sys.exit(1)

    ssearch_big_table = sys.argv[1]
    qstarts_path = sys.argv[2]

    schema = {
        "qseqid": pl.Utf8,
        "sseqid": pl.Utf8,
        "qstart_mapped": pl.Int32,
        "qstart_nuc": pl.Int32,
        "qstart_difference": pl.Int32,
        "qstart": pl.Int32,
        "category": pl.Utf8,
        "recovered_elongation": pl.Float32,
        "raw_recovered_elongation": pl.Int32,
        "qstart_prot": pl.Int32
    }

    columns = [
        "qseqid",
        "sseqid",
        "qstart_mapped",
        "qstart_nuc",
        "qstart_difference",
        "qstart",
        "category",
        "recovered_elongation",
        "raw_recovered_elongation",
        "qstart_prot"
    ]

    qseqid = os.path.basename(ssearch_big_table).split("_big")[0]

    # A query whose subjects produced no alignment at all leaves a header-only
    # file. read_csv raises NoDataError on those, so catch it before parsing and
    # emit a well-formed empty table instead: an empty result is a legitimate
    # outcome, not a failure, and downstream collectFile still needs the file.
    if (os.path.getsize(ssearch_big_table) == 0
            or sum(1 for _ in open(ssearch_big_table)) <= 1):
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
        .select(["qseqid", "sseqid", "qstart_nuc"])
    )

    if table.height == 0:
        pl.DataFrame(schema=schema).write_csv(
            f"{qseqid}_final_res.tsv", separator="\t", include_header=True)
        sys.exit(0)

    qstarts_mapped = pl.read_csv(qstarts_path, separator="\t", has_header=True)

    # ssearch only reports subjects that produced a nucleotide alignment, so its
    # subjects are a subset of qstarts (which lists every subject). Only the
    # reverse — ssearch carrying a subject absent from qstarts — signals a real
    # plumbing bug.
    if table.select(pl.col("sseqid")).n_unique() > qstarts_mapped.select(pl.col("sseqid")).n_unique():
        raise ValueError(
            "The ssearch output has more unique subjects than the qstarts_mapped file "
            "(ssearch subjects should be a subset of qstarts subjects)"
        )

    (
        table
        .join(qstarts_mapped, on="sseqid", how="inner")
        .with_columns(
            qstart_difference=pl.col("qstart_mapped") - pl.col("qstart_nuc"),
            qstart_prot=pl.col("qstart"),
            # qstart == 1 means there is no N-ter extension to recover, so the
            # denominator is 0. Guard it: polars would silently yield inf.
            recovered_elongation=pl.when(pl.col("qstart") > 1)
                .then((((pl.col("qstart") - 1) * 3 + 1) - pl.col("qstart_nuc"))
                      / ((pl.col("qstart") - 1) * 3))
                .otherwise(None),
            raw_recovered_elongation=(((pl.col("qstart") - 1) * 3 + 1) - pl.col("qstart_nuc")),
            category=pl.lit("big")
        )
        .select(columns)
        .cast(schema)
        .write_csv(f"{qseqid}_final_res.tsv", separator="\t", include_header=True)
    )


if __name__ == "__main__":
    main()