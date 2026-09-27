#!/usr/bin/env python3
"""
analyze_small.py - Analyze small N-terminal extension results.

This script classifies annotations based on various criteria like ATG presence,
stop codons, and recovered elongation.

Usage:
    python3 analyze_small.py <small_nter.tsv>

Output:
    Annotated TSV with elongation_status column
"""

import sys
import polars as pl


def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <small_nter.tsv>", file=sys.stderr)
        sys.exit(1)

    small_nter_path = sys.argv[1]

    small = pl.read_csv(small_nter_path, separator="\t", has_header=True)

    atg_not_on_elongated = (
        (pl.col("atg_on_elongated_subject_facing_query_start") == False) |
        ((pl.col("atg_on_elongated_subject_facing_query_start") == True) & (pl.col("stop_inframe") == True))
    )
    
    meth_not_facing_meth = (
        (pl.col("sstart_one") == False) |
        ((pl.col("sstart_one") == True) & (pl.col("meth_on_query_facing_subject_start") == False))
    )

    good_annotation = atg_not_on_elongated & meth_not_facing_meth

    signal_found = (pl.col("sstart_nuc_gt_elongation") == False)
    recovered = (pl.col("raw_recovered_elongation") > 30)

    perfect = good_annotation & signal_found & recovered
    dubious = ~good_annotation & signal_found & recovered
    not_found = ~recovered

    small = (
        small
        .with_columns(
            pl.when(perfect)
            .then(pl.lit("perfect"))
            .when(dubious)
            .then(pl.lit("dubious"))
            .when(not_found)
            .then(pl.lit("not_found"))
            .otherwise(pl.lit("unknown"))
            .alias("elongation_status")
        )
    )

    output_name = small_nter_path.replace(".tsv", "_analyzed.tsv")
    small.write_csv(output_name, separator="\t", include_header=True)


if __name__ == "__main__":
    main()
