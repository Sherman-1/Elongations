#!/usr/bin/env python3
"""
parse_blastp_nter.py - Parse BLASTp results for N-terminal extension analysis.

This script filters BLASTp results to identify homologs, categorizes them
based on sstart position (big vs small), and outputs Parquet files per query.

Usage:
    python3 parse_blastp_nter.py <blast_output.tsv> <coding_coords.tsv> <neighbors_gff_dir>

Output:
    *_small.parquet - Homologs with sstart < 5
    *_big.parquet - Homologs with sstart >= 5
    bad_species.tsv - Species where N-ter is already covered
    full.blast - Filtered BLAST results
    mRNA_species.tsv - mRNA to species mapping
"""

import sys
import glob
import os
import polars as pl
import gff3_parser as gff3


def main():
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <blast_output> <coding_coords> <neighbors_gff_dir>", file=sys.stderr)
        sys.exit(1)

    blast_output = sys.argv[1]
    coding_coords_path = sys.argv[2]
    neighbors_gff_dir = sys.argv[3]

    log = open("log.txt", "w")

    # Iterate through all GFFs in neighbors folder
    mRNA_species = pl.DataFrame({"ID": [], "species": []}, schema={"ID": pl.Utf8, "species": pl.Utf8})

    files = glob.glob(os.path.join(neighbors_gff_dir, "*.gff"))
    log.write(f"Found {len(files)} files\n")

    # Write a TSV file where each mRNA (i.e protein) is associated with its species
    for file in files:
        baseName = os.path.basename(file).split(".")[0]
        log.write(f"Processing {baseName}\n")
        log.write(f"Reading {file}\n")
        mRNA_species.extend(
            (
                pl.from_pandas(gff3.parse_gff3(file, verbose=False, parse_attributes=True))
                .filter(pl.col("Type") == "mRNA")
                .select(pl.col("ID"))
                .cast({"ID": pl.Utf8})
                .with_columns(pl.lit(baseName).alias("species"))
            )
        )

    mRNA_species.write_csv("mRNA_species.tsv", separator="\t", include_header=True)

    # Upstream is the distance between the start of the protein and the start of the gene
    coding_coord = pl.read_csv(coding_coords_path, separator="\t", has_header=True).select(
        ["coding_RNA", "upstream", "Strand"]
    )

    log.write(f"Shape of coding_coord : {coding_coord.shape}\n")

    df = pl.read_csv(blast_output, separator="\t", has_header=True)

    nb_unique_query = len(df.select(["qseqid"]).unique().to_series().to_list())
    log.write(f"Number of unique query before homology filter : {nb_unique_query}\n")
    log.write("We should expect the same number of queries as outputed by Diamond filtered blastp\n")

    df = (
        df
        .with_columns(
            scovhsp=pl.when(pl.col("slen") > 0)
                      .then(((pl.col("send") - pl.col("sstart")).abs() + 1)
                            / pl.col("slen") * 100)
                      .otherwise(0)
        )
        .filter(
            (pl.col("evalue") <= 1e-5)
            & (pl.col("qcovhsp") >= 70)
            & (pl.col("ppos") > 50)
        )
        .join(coding_coord, left_on="sseqid", right_on="coding_RNA", how="inner")
        .join(mRNA_species, left_on="sseqid", right_on="ID", how="inner")
    )

    df.write_csv("full.blast", separator="\t", include_header=True)
    log.write(f"Blasts results filtered for homology 'full.blast': {df.shape}\n")

    nb_unique_query = len(df.select(["qseqid"]).unique().to_series().to_list())
    log.write(f"Number of unique query after homology filter : {nb_unique_query}\n")

    # Search every query and its associated taxid that has a hit with
    # qstart < 20. This step filters out species where the nter of the query is
    # seen in at least one protein.
    bad_seqs = df.group_by("qseqid").agg(
        [pl.col("species").filter(pl.col("qstart") < 20).unique().alias("bad_species")]
    ).select(["qseqid", "bad_species"])

    bad_seqs.with_columns(pl.col("bad_species").list.join(separator=";")).write_csv(
        "bad_species.tsv", separator="\t", include_header=True
    )

    # Then, we filter out the bad species from the full blast result
    # To only keep the species where the nter is not seen in any protein
    filtered = (
        df.join(bad_seqs, on="qseqid", how="inner")
        .filter(pl.col("species").is_in(pl.col("bad_species")).not_())
        .filter(pl.col("scovhsp") >= 70)
        .select(pl.exclude(["bad_species"]))
    )

    # Supposedly, there is no more row with a qstart < 20
    # Test it, and raise an error if it is not the case
    if filtered.filter(pl.col("qstart") < 20).shape[0] > 0:
        filtered.write_csv("error.blast", separator="\t", include_header=True)
        raise ValueError("There are still some proteins with qstart < 20")

    filtered.write_csv("filtered.blast", separator="\t", include_header=True)
    log.write(f"Blasts results filtered again on bad species 'filtered.blast': {filtered.shape}\n")

    nb_unique_query = len(filtered.select(["qseqid"]).unique().to_series().to_list())
    log.write(f"Number of unique query after bad species filter : {nb_unique_query}\n")

    big = filtered.filter((pl.col("sstart") >= 5)).select(["qseqid", "sseqid", "qstart", "sstart"])

    small = filtered.filter((pl.col("sstart") < 5)).select(
        ["qseqid", "sseqid", "qstart", "sstart", "upstream", "Strand", "slen"]
    )

    # In theory, big + small = filtered
    if filtered.shape[0] != big.shape[0] + small.shape[0]:
        log.write(f"Shape of filtered : {filtered.shape}\n")
        log.write(f"Shape of big : {big.shape}\n")
        log.write(f"Shape of small : {small.shape}\n")
        raise ValueError("big + small != filtered")

    nb_unique_query_big = big.select(["qseqid"]).unique().to_series().to_list()
    nb_unique_query_small = small.select(["qseqid"]).unique().to_series().to_list()

    log.write(f"Number of unique query with big sstart : {len(nb_unique_query_big)}\n")
    log.write(f"Number of unique query with small sstart : {len(nb_unique_query_small)}\n")

    # Are there some overlap between big and small?
    inter = len(set(nb_unique_query_big).intersection(set(nb_unique_query_small)))
    log.write(f"Overlap between big and small : {inter}\n")

    # Does nb_unique_query = nb_unique_query_big + nb_unique_query_small ?
    log.write(
        f"Does nb_unique_query = nb_unique_query_big + nb_unique_query_small ? "
        f"{nb_unique_query == len(nb_unique_query_big) + len(nb_unique_query_small)}\n"
    )
    log.write(f"{nb_unique_query} = {len(nb_unique_query_big)} + {len(nb_unique_query_small)}\n")

    big.write_csv("big.blast", separator="\t", include_header=True)
    small.write_csv("small.blast", separator="\t", include_header=True)

    for (query,), big_group in big.group_by("qseqid"):
        big_group.write_parquet(f"{query}_big.parquet", compression="snappy")

    for (query,), small_group in small.group_by("qseqid"):
        small_group.write_parquet(f"{query}_small.parquet", compression="snappy")

    log.write(f"Shape of full after writings : {df.shape}")
    log.close()


if __name__ == "__main__":
    main()
