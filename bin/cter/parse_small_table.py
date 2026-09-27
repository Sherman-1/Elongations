#!/usr/bin/env python3
"""
parse_small_table.py - Parse ssearch results for small send C-terminal candidates.

This script performs analysis of nucleotide alignments to detect
potential C-terminal extensions.

Usage:
    python3 parse_small_table.py <ssearch_elong.tsv> <ssearch_elong.aln> <ssearch_short.tsv> \
                                  <tfastx.tsv> <thresholds.tsv> <elongated_subjects.faa>

Output:
    {qseqid}_final_res.tsv
"""

import sys
import os
import csv
import re
from collections import defaultdict
from Bio import SeqIO
import polars as pl


# =============================
#         Functions
# =============================

def read_thresholds(thresholds_path: str, log):
    """
    Reads tabular file, schema is as follow:
    "subject elongated_length elongation slend_mapped_nuc qend qlend slen"
    """
    thresholds = defaultdict(dict)
    with open(thresholds_path, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            log.write(f"Row : {row}\n")

            if row == []:
                continue

            subject = row[0]
            elongated_length = int(row[1])
            elongation = int(row[2])
            slend_mapped_nuc = int(row[3])
            qend = int(row[4])
            qlend = int(row[5])
            slen_prot = int(row[6])

            thresholds[subject] = {
                "elongated_length": elongated_length,
                "elongation": elongation,
                "slend_mapped_nuc": slend_mapped_nuc,
                "qend_blastp": qend,
                "qlend_blastp": qlend,
                "slen_prot": slen_prot
            }

    return thresholds


def parse_ssearch_output(ssearch_aln_path: str, side: str, log):
    """
    This function parses the '-m A' format output of the ssearch36 program
    and returns a dict of sequence + informations, one per subject sequence.
    """
    assert side in ["nter", "cter"]

    with open(ssearch_aln_path, "r") as f:
        lines = f.readlines()

    flag = False
    skip_count = -1
    hsp_number = -1
    subject = None

    pattern_overlap = re.compile(r'in (\d+) nt overlap')
    pattern_bitscore = re.compile(r'bits:\s+([0-9]+\.[0-9]+)')

    datas = defaultdict(dict)
    records = defaultdict(dict)
    qstarts = defaultdict(dict)
    qends = defaultdict(dict)
    overlaps = defaultdict(dict)
    bitscores = defaultdict(dict)

    for line in lines:
        line = line.strip()

        if not line:
            flag = False
            continue

        if skip_count > 0:
            if skip_count == 3:
                bitscore = float(pattern_bitscore.search(line).group(1))
                bitscores[subject][hsp_number] = bitscore
                log.write(f"Bitscore for {subject} hsp {hsp_number}: {bitscore}\n")
            elif skip_count == 2:
                overlaps[subject][hsp_number] = int(pattern_overlap.search(line).group(1))
            skip_count -= 1
            continue

        if line.startswith(">>") or line == ">--":
            flag = True
            if line.startswith(">>"):
                subject = line.split()[0][2:]
                hsp_number = max(datas[subject].keys()) + 1 if subject in datas else 0
            else:
                hsp_number += 1

            skip_count = 3

            if subject not in datas:
                datas[subject] = {}

            datas[subject].update({hsp_number: {"positions": [], "subject": [], "query": []}})
            qstarts[subject].update({hsp_number: 999999})
            qends[subject].update({hsp_number: 0})
            continue

        if flag:
            line_parts = line.split()
            subject_nucleotide_number = int(line_parts[1])
            query_nucleotide = line_parts[2]
            subject_nucleotide = line_parts[3]

            datas[subject][hsp_number]["positions"].append(subject_nucleotide_number)
            datas[subject][hsp_number]["subject"].append(subject_nucleotide)
            datas[subject][hsp_number]["query"].append(query_nucleotide)

            qstart = int(line_parts[0])
            qstarts[subject][hsp_number] = min(qstarts[subject][hsp_number], qstart)
            qends[subject][hsp_number] = max(qends[subject][hsp_number], qstart)

    for subject in datas:
        hsp_number, _ = max(bitscores[subject].items(), key=lambda x: x[1])

        if side == "cter":
            qend = qends[subject][hsp_number]
            records[subject] = {
                "subject": datas[subject][hsp_number]["subject"],
                "query": datas[subject][hsp_number]["query"],
                "positions": datas[subject][hsp_number]["positions"],
                "qend": qend
            }

    log.write(f"Records : {records}\n")
    return dict(records)


def main():
    if len(sys.argv) != 7:
        print(f"Usage: {sys.argv[0]} <ssearch_elong.tsv> <ssearch_elong.aln> <ssearch_short.tsv> "
              f"<tfastx.tsv> <thresholds.tsv> <elongated_subjects.faa>", file=sys.stderr)
        sys.exit(1)

    ssearch_elong_table = sys.argv[1]
    ssearch_elong_aln = sys.argv[2]
    ssearch_short_table = sys.argv[3]
    tfastx_table = sys.argv[4]
    thresholds_file = sys.argv[5]
    elongated_subjects_faa = sys.argv[6]

    log = open("log.txt", "w")
    output = open("output.txt", "w")

    # Check if files are not empty
    if (os.stat(ssearch_elong_table).st_size == 0 or
            os.stat(ssearch_elong_aln).st_size == 0 or
            os.stat(thresholds_file).st_size == 0 or
            os.stat(tfastx_table).st_size == 0):
        log.write("One of the input files is empty\n")
        log.close()
        output.close()
        return

    qseqid = os.path.basename(ssearch_elong_table).split("_complete")[0]

    elongated_proteins = SeqIO.to_dict(SeqIO.parse(elongated_subjects_faa, "fasta"))

    thresholds = read_thresholds(thresholds_file, log)

    thresholds_df = pl.DataFrame({
        "sseqid": list(thresholds.keys()),
        "elongated_length": [v["elongated_length"] for v in thresholds.values()],
        "elongation": [v["elongation"] for v in thresholds.values()],
        "slend_mapped_nuc": [v["slend_mapped_nuc"] for v in thresholds.values()],
        "qend_blastp": [v["qend_blastp"] for v in thresholds.values()],
        "qlend_blastp": [v["qlend_blastp"] for v in thresholds.values()],
        "slen_prot": [v["slen_prot"] for v in thresholds.values()]
    })

    # Read alignment tables
    table_nuc_elong = (
        pl.read_csv(ssearch_elong_table, separator="\t", has_header=True,
                    schema_overrides={"evalue": pl.Float32})
        .sort(["qseqid", "sseqid", "bitscore"], descending=[False, False, True])
        .unique(subset=["qseqid", "sseqid"], keep="first")
        .select(["qseqid", "sseqid", "send_nuc_elong", "qend_nuc_elong"])
    )

    if table_nuc_elong.select(pl.col("sseqid")).n_unique() != table_nuc_elong.shape[0]:
        raise ValueError("There are duplicates in the table_nuc_elong")

    table_nuc_short = (
        pl.read_csv(ssearch_short_table, separator="\t", has_header=True,
                    schema_overrides={"evalue": pl.Float32})
        .sort(["qseqid", "sseqid", "bitscore"], descending=[False, False, True])
        .unique(subset=["qseqid", "sseqid"], keep="first")
        .with_columns(
            qlend_nuc_short=(pl.col("qlen_nuc_short") - pl.col("qend_nuc_short"))
        )
        .select(["sseqid", "qlend_nuc_short", "send_nuc_short", "qend_nuc_short"])
    )

    if table_nuc_short.select(pl.col("sseqid")).n_unique() != table_nuc_short.shape[0]:
        raise ValueError("There are duplicates in the table_nuc_short")

    table_tfastx = (
        pl.read_csv(tfastx_table, separator="\t", has_header=True,
                    schema_overrides={"evalue": pl.Float32})
        .sort(["qseqid", "sseqid", "bitscore"], descending=[False, False, True])
        .unique(subset=["qseqid", "sseqid"], keep="first")
        .with_columns(
            qlend_tfastx=(pl.col("qlen_tfastx") - pl.col("qend_tfastx"))
        )
        .select(["sseqid", "qlend_tfastx", "send_tfastx", "qend_tfastx"])
    )

    if table_tfastx.select(pl.col("sseqid")).n_unique() != table_tfastx.shape[0]:
        raise ValueError("There are duplicates in the table_tfastx")

    diff = (
        table_tfastx
        .join(table_nuc_elong, on=["sseqid"], how="inner")
        .join(table_nuc_short, on=["sseqid"], how="inner")
        .join(thresholds_df, on=["sseqid"], how="inner")
        .with_columns(
            qend_prot_diff=(pl.col("qend_tfastx") - pl.col("qend_blastp")),
            qend_nuc_diff=(pl.col("qend_nuc_elong") - pl.col("qend_nuc_short"))
        )
    ).select(["qseqid", "sseqid", "qlend_tfastx", "qlend_blastp", "qend_prot_diff",
              "qend_nuc_diff", "qlend_nuc_short", "qend_nuc_elong", "qend_nuc_short"])

    diff.write_csv("diff.tsv", separator="\t", include_header=True)

    log.write("Reading alignments ...\n")
    alignments = parse_ssearch_output(ssearch_elong_aln, "cter", log)

    log.write("Checking for subjects in the table, thresholds and alignments ...\n")
    table_subjects = table_nuc_elong.select("sseqid").unique().to_series().to_list()
    threshold_subjects = thresholds_df.select("sseqid").unique().to_series().to_list()
    alignment_subjects = list(alignments.keys())

    if not (set(table_subjects) == set(threshold_subjects) == set(alignment_subjects)):
        log.write("Error, the subjects in the table, thresholds and alignments are not the same\n")
        uniques = set(table_subjects) ^ set(threshold_subjects) ^ set(alignment_subjects)
        log.write(f"Difference : {uniques}\n")

    try:
        workdir = os.getcwd().split("work/")[1]
    except IndexError:
        workdir = os.getcwd()

    columns = [
        "qseqid","sseqid", "qend_nuc_elong", "qend_nuc_short", "qlend_nuc_short",
        "qlend_tfastx", "qlend_blastp", "qend_prot_diff", "qend_nuc_diff",
        "workdir", "recovered_elongation", "raw_recovered_elongation", "category"
    ]

    schema = {
        "qseqid": pl.Utf8,
        "sseqid": pl.Utf8,
        "qend_nuc_elong": pl.Int32,
        "qend_nuc_short": pl.Int32,
        "qlend_nuc_short": pl.Int32,
        "qlend_tfastx": pl.Int32,
        "qlend_blastp": pl.Int32,
        "qend_prot_diff": pl.Int32,
        "qend_nuc_diff": pl.Int32,
        "workdir": pl.Utf8,
        "recovered_elongation": pl.Float32,
        "raw_recovered_elongation": pl.Int32,
        "category": pl.Utf8
    }

    (
        diff
        .join(thresholds_df, on=["sseqid"], how="inner")
        .with_columns(
            workdir=pl.lit(workdir),
            recovered_elongation=(
                (pl.col("qend_nuc_elong") - ((pl.col("qend_blastp") - 1) * 3 + 1))
                / (pl.col("qlend_blastp") * 3)
            ),
            raw_recovered_elongation=(
                pl.col("qend_nuc_elong") - ((pl.col("qend_blastp") - 1) * 3 + 1)
            ),
            category=pl.lit("small")
        )
        .select(columns)
        .cast(schema)
    ).write_csv(f"{qseqid}_final_res.tsv", separator="\t", include_header=True)

    log.close()
    output.close()


if __name__ == "__main__":
    main()
