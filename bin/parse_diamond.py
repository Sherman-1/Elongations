#!/usr/bin/env python3
"""
parse_diamond.py - Parse Diamond BLAST output to identify N-ter and C-ter candidates.

This script analyzes Diamond output against NR database to identify proteins
where less than 5% of eukaryotic homologs have alignments covering the
N-terminal or C-terminal regions.

Usage:
    python3 parse_diamond.py <diamond_output.tsv> <strain2species.csv> <eukaryotes.csv>

Output:
    nter_candidates_<num>.tsv - Proteins with potential N-terminal extensions
    cter_candidates_<num>.tsv - Proteins with potential C-terminal extensions
    statistics_<num>.tsv - Full statistics per protein
"""

import sys
import os
import csv
import polars as pl


def main():
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <diamond_df> <strain2species> <eukaryotes>", file=sys.stderr)
        sys.exit(1)

    df_path = sys.argv[1]
    strain2species_path = sys.argv[2]
    eukaryotes_path = sys.argv[3]

    bn = os.path.basename(df_path)
    number = bn.split(".")[0]

    # Load eukaryotes taxid mapping
    with open(eukaryotes_path) as f:
        reader = csv.reader(f, delimiter=";")
        eukaryotes = {int(row[0]): int(row[1]) for row in reader if row[0] != "" and row[1] != ""}

    # Load strain to species mapping
    with open(strain2species_path) as f:
        reader = csv.reader(f, delimiter=";")
        species_map = {int(row[0]): int(row[1]) for row in reader if row[0] != "" and row[1] != ""}

    lazy = (
        pl.scan_csv(
            df_path,
            separator="\t",
            has_header=False,
            new_columns=["qseqid", "sseqid", "qlen", "qstart", "qlend", "staxids"],
            schema_overrides={
                "qseqid": pl.Utf8,
                "sseqid": pl.Utf8,
                "qlen": pl.Int64,
                "qstart": pl.Int32,
                "qlend": pl.Int64,
                "staxids": pl.Utf8
            }
        )
        .with_columns(pl.col("staxids").str.split(";"))
        .explode("staxids")
        .with_columns(pl.col("staxids").cast(pl.Int32))
        .with_columns([
            pl.col("staxids").replace(eukaryotes, default=0).alias("euk"),
            pl.col("staxids").replace(species_map, default=0).alias("species")
        ])
        .filter(pl.col("euk") == 1)
        .filter(pl.col("species") != 0)
        .group_by("qseqid").agg([
            pl.col("species").n_unique().alias("nb_unique_species").cast(pl.Int32),
            pl.col("species").filter(pl.col("qstart") <= 15).n_unique().alias("nter_q_15").cast(pl.Int32),
            pl.col("species").filter(pl.col("qlend") <= 15).n_unique().alias("cter_q_15").cast(pl.Int32),
            pl.col("qlen").first().alias("qlen"),
            pl.col("qstart").median().alias("qstart_median"),
            pl.col("qstart").std().alias("qstart_std"),
            pl.col("qlend").median().alias("cter_median"),
            pl.col("qlend").std().alias("cter_std")
        ])
        .with_columns([
            (pl.col("nter_q_15") / pl.col("nb_unique_species") * 100).alias("percent_nter_q_15").cast(pl.Float32),
            (pl.col("cter_q_15") / pl.col("nb_unique_species") * 100).alias("percent_cter_q_15").cast(pl.Float32),
        ])
    )

    df = lazy.collect()
    df.write_csv(f"statistics_{number}.tsv", separator="\t")
    df.filter(pl.col("percent_nter_q_15") <= 5).select("qseqid").write_csv(
        f"nter_candidates_{number}.tsv", separator="\t", include_header=False
    )
    df.filter(pl.col("percent_cter_q_15") <= 5).select("qseqid").write_csv(
        f"cter_candidates_{number}.tsv", separator="\t", include_header=False
    )


if __name__ == "__main__":
    main()
