#!/usr/bin/env python3
"""
check_annotation.py - Check for duplicate mRNA IDs across GFF files.

This script reads all GFF files in the input directories and checks
if any mRNA ID appears in multiple files, which would indicate
annotation conflicts.

Usage:
    python3 check_annotation.py <gff_pattern>

Example:
    python3 check_annotation.py "input/*/*.gff"

Output:
    Raises exception if duplicates found, otherwise exits silently.
"""

import sys
import glob
import polars as pl
import gff3_parser as gff3


def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <gff_glob_pattern>", file=sys.stderr)
        sys.exit(1)

    pattern = sys.argv[1]
    files = glob.glob(pattern)

    lst = []
    print("Reading files")
    for file in files:
        lst.append(
            pl.from_pandas(gff3.parse_gff3(file, verbose=False, parse_attributes=True))
            .filter(
                (pl.col("Type") == "mRNA")
            )
            .select([
                "ID"
            ])
            .to_series()
            .to_list()
        )

    print("Files read")

    flag = False

    for i, sublist1 in enumerate(lst, start=1):
        for j, sublist2 in enumerate(lst, start=1):
            if i != j:  # Skip comparing the same list
                intersection = set(sublist1).intersection(sublist2)
                if intersection:
                    print(f"Lists {i} and {j} have common items: {intersection}")
                    flag = True

    if flag:
        raise Exception("Common mRNAs identified in different files. Please check the input files.")


if __name__ == "__main__":
    main()
