#!/usr/bin/env python3
"""
get_gene_coords.py - Calculate upstream/downstream coordinates from GFF annotations.

This script reads a GFF3 file and calculates the distance between CDS and gene
boundaries for each mRNA, outputting a TSV with upstream/downstream distances
adjusted for strand orientation.

Usage:
    python3 get_gene_coords.py <gff_file>

Output:
    gene_coords.tsv
"""

import sys
import polars as pl
import gff3_parser as gff3


def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <gff_file>", file=sys.stderr)
        sys.exit(1)

    gff_path = sys.argv[1]

    gff = pl.from_pandas(gff3.parse_gff3(gff_path, verbose=False, parse_attributes=True))

    (
        gff
        .filter(
            (pl.col("Type") == "CDS")
        )
        .group_by("Parent").agg([
            pl.col("Start").min().alias("Start"),
            pl.col("End").max().alias("End"),
            pl.col("Strand").first().alias("Strand")
        ])
        .rename({
            "Parent": "coding_RNA",
            "Start": "protein_start",
            "End": "protein_end"
        })
        .join(
            other=(
                gff
                .filter(
                    (pl.col("Type") == "mRNA")
                )
                .select([
                    "Parent",
                    "ID",
                ])
                .join(
                    gff.filter(pl.col("Type") == "gene").select(["ID", "Parent", "Start", "End"]),
                    left_on="Parent",
                    right_on="ID",
                    suffix="_gene",
                    how="inner"
                )
                .drop("Parent_gene")  # Drop gene parent, i.e nobody
                .rename({
                    "Parent": "gene_ID",
                    "ID": "mRNA_ID",
                    "Start": "gene_start",
                    "End": "gene_end"
                })
            ),
            left_on="coding_RNA",
            right_on="mRNA_ID",
            how="inner",
            suffix="_mRNA"
        )
        .cast({
            "protein_start": pl.Int32,
            "protein_end": pl.Int32,
            "gene_start": pl.Int32,
            "gene_end": pl.Int32
        })
        .with_columns([
            (pl.col("protein_start") - pl.col("gene_start")).alias("upstream"),
            (pl.col("gene_end") - pl.col("protein_end")).alias("downstream")
        ])
        .with_columns(
            pl.when(pl.col("Strand") == "-")
            .then(
                pl.struct([
                    pl.col("downstream").alias("upstream"),
                    pl.col("upstream").alias("downstream")
                ])
            )
            .otherwise(
                pl.struct([
                    pl.col("downstream").alias("downstream"),
                    pl.col("upstream").alias("upstream")
                ])
            )
            .alias("new_coords")
        )
        .drop(["downstream", "upstream"])
        .unnest("new_coords")
    ).write_csv("gene_coords.tsv", separator="\t", include_header=True)


if __name__ == "__main__":
    main()
