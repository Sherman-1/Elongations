#!/usr/bin/env python3
"""
parse_small_table.py - Parse ssearch results for small sstart N-terminal candidates.

This script performs complex analysis of nucleotide alignments to detect
potential N-terminal extensions, including:
- Gap analysis in alignments
- Start/stop codon detection (in-frame)
- Recovery of elongation calculation

Usage:
    python3 parse_small_table.py <ssearch_elong.tsv> <ssearch_elong.aln> <ssearch_short.tsv> \
                                  <tfastx.tsv> <thresholds.tsv> <elongated_subjects.fna> <query.faa>

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

def read_thresholds(thresholds_path: str):
    """
    Reads tabular file, schema is as follow:
    subject threshold elongation sstart_mapped_nuc qstart elongated_length sstart
    The file is headerless
    """
    thresholds = defaultdict(dict)
    with open(thresholds_path, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if row == []:
                continue

            subject = row[0]
            threshold = int(row[1])
            elongation = int(row[2])
            sstart_mapped_nuc = int(row[3])
            qstart_blastp = int(row[4])
            elongated_subject_length = int(row[5])
            sstart_blastp = int(row[6])

            thresholds[subject]["threshold"] = threshold
            thresholds[subject]["elongation"] = elongation
            thresholds[subject]["sstart_mapped_nuc"] = sstart_mapped_nuc
            thresholds[subject]["qstart_blastp"] = qstart_blastp
            thresholds[subject]["elongated_length"] = elongated_subject_length
            thresholds[subject]["sstart_blastp"] = sstart_blastp

    return thresholds


def parse_ssearch_output(ssearch_aln_path: str, side: str, log):
    """
    This function parses the '-a -mA' format output of the ssearch36 program
    and returns a dict of sequence + informations, one per subject sequence.
    If multiple HSPs are present for a given query-subject pair, only the one
    with the highest bitscore is kept.
    """
    assert side in ["nter", "cter"]

    with open(ssearch_aln_path, "r") as f:
        lines = f.readlines()

    flag = False
    skip_count = -1
    hsp_number = -1
    subject = None

    #pattern_headers = re.compile(r'^>\w+')
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

        if side == "nter":
            qstart = qstarts[subject][hsp_number]
            records[subject] = {
                "subject": datas[subject][hsp_number]["subject"],
                "query": datas[subject][hsp_number]["query"],
                "positions": datas[subject][hsp_number]["positions"],
                "qstart": qstart
            }

            for hsp in qstarts[subject]:
                log.write(f"Subject : {subject}, hsp : {hsp}, qstart : {qstarts[subject][hsp]}, bitscore : {bitscores[subject][hsp]}\n")

    log.write(f"Records : {records}\n")
    return dict(records)


def check_atg_gaps(alignments: dict, thresholds: dict, subject: str, side: str, log):
    """
    Count gaps between the end (or start) of the sequence up to the threshold position.
    """
    assert side in ["nter", "cter"]

    positions = alignments[subject]["positions"]
    query = alignments[subject]["query"]
    sequence = alignments[subject]["subject"]

    threshold = thresholds[subject]["threshold"]
    log.write(f"Threshold for {subject} : {threshold}\n")

    gaps_subject = 0
    gaps_query = 0

    if side == "nter":
        for position, subj_nucleotide, query_nucleotide in zip(positions, sequence, query):
            if position > threshold:
                break
            if subj_nucleotide == "-":
                gaps_subject += 1
                continue
            if query_nucleotide == "-":
                gaps_query += 1
                continue

        return {
            "sseqid": subject,
            "gaps_query": gaps_query,
            "gaps_subject": gaps_subject,
            "qstart_nuc_elong": alignments[subject]["qstart"],
            "sstart_nuc_gt_elongation": False
        }


def check_inframe_codon(sequence, start, stop, codon_set):
    """
    Check for the presence of specific codons in the sequence between start and stop positions.
    """
    assert type(sequence) == str
    assert type(start) == int
    assert type(stop) == int

    if start <= stop:
        return False, -1

    for i in range(start, stop, -3):
        codon = sequence[i - 3:i].upper()
        if codon in codon_set:
            return True, i
    return False, -1


def main():
    if len(sys.argv) != 8:
        print(f"Usage: {sys.argv[0]} <ssearch_elong.tsv> <ssearch_elong.aln> <ssearch_short.tsv> "
              f"<tfastx.tsv> <thresholds.tsv> <elongated_subjects.fna> <query.faa>", file=sys.stderr)
        sys.exit(1)

    ssearch_elong_table = sys.argv[1]
    ssearch_elong_aln = sys.argv[2]
    ssearch_short_table = sys.argv[3]
    tfastx_table = sys.argv[4]
    thresholds_file = sys.argv[5]
    elongated_subjects_fna = sys.argv[6]
    query_faa = sys.argv[7]

    log = open("log.txt", "w")
    output = open("output.txt", "w")

    if (os.stat(ssearch_elong_table).st_size == 0 or
            os.stat(ssearch_elong_aln).st_size == 0 or
            os.stat(thresholds_file).st_size == 0 or
            os.stat(tfastx_table).st_size == 0):
        log.write("One of the input files is empty\n")
        log.close()
        output.close()
        return

    qseqid = os.path.basename(ssearch_elong_table).split("_complete")[0]

    # Read subject sequences
    elongated_subjects_sequences = SeqIO.to_dict(SeqIO.parse(elongated_subjects_fna, "fasta"))

    # Read query sequence
    query = SeqIO.read(query_faa, "fasta")
    query = str(query.seq)

    # Extract thresholds for each subject sequence
    thresholds = read_thresholds(thresholds_file)

    thresholds_df = pl.DataFrame({
        "sseqid": list(thresholds.keys()),
        "threshold": [v["threshold"] for v in thresholds.values()],
        "elongation": [v["elongation"] for v in thresholds.values()],
        "sstart_mapped_nuc": [v["sstart_mapped_nuc"] for v in thresholds.values()],
        "qstart_blastp": [v["qstart_blastp"] for v in thresholds.values()],
        "elongated_length": [v["elongated_length"] for v in thresholds.values()],
        "sstart_blastp": [v["sstart_blastp"] for v in thresholds.values()]
    })

    # Read alignment tables
    table_nuc_elong = (
        pl.read_csv(ssearch_elong_table, separator="\t", has_header=True,
                    schema_overrides={"evalue": pl.Float32})
        .sort(["qseqid", "sseqid", "bitscore"], descending=[False, False, True])
        .unique(subset=["qseqid", "sseqid"], keep="first")
        .select(["qseqid", "sseqid", "qstart_nuc_elong", "sstart_nuc_elong"])
    )

    if table_nuc_elong.select(pl.col("sseqid")).n_unique() != table_nuc_elong.shape[0]:
        raise ValueError("There are duplicates in the table_nuc_elong")

    table_nuc_short = (
        pl.read_csv(ssearch_short_table, separator="\t", has_header=True,
                    schema_overrides={"evalue": pl.Float32})
        .sort(["qseqid", "sseqid", "bitscore"], descending=[False, False, True])
        .unique(subset=["qseqid", "sseqid"], keep="first")
        .select(["sseqid", "qstart_nuc_short", "sstart_nuc_short"])
    )

    if table_nuc_short["sseqid"].n_unique() != table_nuc_short.shape[0]:
        raise ValueError("There are duplicates in the table_nuc_short")

    table_tfastx = (
        pl.read_csv(tfastx_table, separator="\t", has_header=True,
                    schema_overrides={"evalue": pl.Float32})
        .sort(["qseqid", "sseqid", "bitscore"], descending=[False, False, True])
        .unique(subset=["qseqid", "sseqid"], keep="first")
        .select(["sseqid", "qstart_tfastx", "sstart_tfastx"])
    )

    if table_tfastx["sseqid"].n_unique() != table_tfastx.shape[0]:
        raise ValueError("There are duplicates in the table_tfastx")

    diff = (
        table_tfastx
        .join(table_nuc_elong, on=["sseqid"], how="inner")
        .join(table_nuc_short, on=["sseqid"], how="inner")
        .join(thresholds_df, on=["sseqid"], how="inner")
        .with_columns(
            qstart_prot_diff=(pl.col("qstart_blastp") - pl.col("qstart_tfastx")),
            qstart_nuc_diff=(pl.col("qstart_nuc_short") - pl.col("qstart_nuc_elong"))
        )
    ).select(["qseqid", "sseqid", "qstart_tfastx", "qstart_blastp", "qstart_prot_diff",
              "qstart_nuc_diff", "qstart_nuc_elong", "qstart_nuc_short"])

    # Read alignments raw
    alignments = parse_ssearch_output(ssearch_elong_aln, "nter", log)

    # Assert good format
    for subject in alignments:
        qstart_nuc_elong = table_nuc_elong.filter(pl.col("sseqid") == subject).select("qstart_nuc_elong").item()
        if qstart_nuc_elong != alignments[subject]["qstart"]:
            log.write(f"Error, qstart_nuc_elong in table and in alignment are not the same for subject {subject}\n")
            raise ValueError("qstart_nuc_elong in table and in alignment are not the same")

    table_subjects = table_nuc_elong.select("sseqid").unique().to_series().to_list()
    threshold_subjects = list(thresholds.keys())
    alignment_subjects = list(alignments.keys())

    if not (set(table_subjects) == set(threshold_subjects) == set(alignment_subjects)):
        log.write("Error, the subjects in the table, thresholds and alignments are not the same\n")
        uniques = set(table_subjects) ^ set(threshold_subjects) ^ set(alignment_subjects)
        log.write(f"Difference : {uniques}\n")
        raise ValueError("Subjects are not the same in table, thresholds and alignments")

    # Compute results
    results = []

    for subject in alignments:
        log.write(f"Processing subject {subject}\n")

        sstart_nuc_elong = table_nuc_elong.filter((pl.col("sseqid") == subject)).select("sstart_nuc_elong").item()
        qstart_nuc_elong = table_nuc_elong.filter((pl.col("sseqid") == subject)).select("qstart_nuc_elong").item()

        qstart_blastp = thresholds[subject]["qstart_blastp"]
        sstart_blastp = thresholds[subject]["sstart_blastp"]

        qstart_tfastx = table_tfastx.filter(pl.col("sseqid") == subject).select("qstart_tfastx").item()
        sstart_tfastx = table_tfastx.filter(pl.col("sseqid") == subject).select("sstart_tfastx").item()

        elongation_length = thresholds[subject]["elongation"]

        # Check methionine in query facing subject start
        if sstart_blastp == 1:
            meth_query_facing_subject_start = query[qstart_blastp - 1] == "M"
            sstart_one = True
        else:
            meth_query_facing_subject_start = True
            sstart_one = False

        assert type(sstart_one) == bool

        # If sstart_nuc_elong is greater than the elongated 
        if sstart_nuc_elong > elongation_length:
            buffer_dict = {
                "sseqid": subject,
                "gaps_query": -1,
                "gaps_subject": -1,
                "qstart_nuc_elong": qstart_nuc_elong,
                "sstart_nuc_gt_elongation": True,
                "stop_inframe": False,
                "stop_inframe_position": -1,
                "start_inframe": True,
                "start_inframe_position": -1,
                "atg_on_elongated_subject_facing_query_start": True,
                "meth_on_query_facing_subject_start": meth_query_facing_subject_start,
                "sstart_one": sstart_one,
                "analysable": False,
                "atg_codons_upstream": -1,
                "atg_explains": None
            }

            log.write(f"Subject {subject} has sstart_sw > threshold\n")
            log.write(f"Buffer dict : {buffer_dict}\n")
            results.append(buffer_dict)

        else:
            buffer_dict = check_atg_gaps(alignments, thresholds, subject, "nter", log)
            assert buffer_dict["sstart_nuc_gt_elongation"] == False
            gaps_query = buffer_dict["gaps_query"]
            #resuls.append(buffer_dict)

            subject_seq = str(elongated_subjects_sequences[subject].seq)

            # Check for start codon on subject facing start on query
            subject_alt_start_position = elongation_length - ((qstart_blastp - 1) * 3 - (sstart_blastp - 1) * 3) - gaps_query + 1
            is_start_on_subject_facing_start_on_query = elongated_subjects_sequences[subject].seq[subject_alt_start_position - 1:subject_alt_start_position + 3 - 1].upper() == "ATG"

            start_position = elongation_length - 2
            stop_position = sstart_nuc_elong - 2

            log.write("Elongated subject sequence near position given by elongation_length - 2\n")
            small_seq = subject_seq[start_position - 10:start_position + 10]
            log.write(f"Small sequence : {small_seq}\n")

            putative_start_codon = subject_seq[start_position + 3 - 1:start_position + 6 - 1].upper()
            log.write(f"Putative start codon : {putative_start_codon}\n")
            assert putative_start_codon == "ATG", (
                f"{subject}: expected ATG at {start_position + 3} in the elongated "
                f"subject, found {putative_start_codon!r} (context: {small_seq!r})")
            log.write(f"Start codon found at position {start_position + 3}\n")

            is_stop_inframe, inframe_stop_position = check_inframe_codon(subject_seq, start_position - 1, stop_position - 1, {"TAA", "TAG", "TGA"})
            is_start_inframe, inframe_start_position = check_inframe_codon(subject_seq, start_position - 1, stop_position - 1, {"ATG"})

            if is_start_inframe:
                atg_codons_upstream = (elongation_length - inframe_start_position) // 3 + 1
                atg_explains = atg_codons_upstream / (qstart_blastp - 1) if qstart_blastp > 1 else None
            else:
                atg_codons_upstream = -1
                atg_explains = None

            assert type(is_stop_inframe) == bool
            assert type(is_start_inframe) == bool
            assert type(meth_query_facing_subject_start) == bool
            assert type(is_start_on_subject_facing_start_on_query) == bool
            assert type(sstart_one) == bool

            buffer_dict.update({
                "meth_on_query_facing_subject_start": meth_query_facing_subject_start,
                "sstart_one": sstart_one,   
                "stop_inframe": is_stop_inframe,
                "stop_inframe_position": inframe_stop_position,
                "start_inframe": is_start_inframe,
                "start_inframe_position": inframe_start_position,
                "atg_on_elongated_subject_facing_query_start": is_start_on_subject_facing_start_on_query,
                "analysable": True,
                "atg_codons_upstream": atg_codons_upstream,
                "atg_explains": atg_explains
            })

            results.append(buffer_dict)

    results = pl.DataFrame(results)
    results.write_csv("results_nter.tsv", separator="\t", include_header=True)
    diff.write_csv("diff.tsv", separator="\t", include_header=True)

    try:
        workdir = os.getcwd().split("work/")[1]
    except IndexError:
        workdir = os.getcwd()

    columns = [
        'qseqid', 'sseqid', 'gaps_query', 'gaps_subject', 'qstart_nuc_elong',
        'sstart_nuc_gt_elongation', 'stop_inframe', 'stop_inframe_position',
        'start_inframe', 'start_inframe_position', 'atg_on_elongated_subject_facing_query_start',
        'meth_on_query_facing_subject_start', 'sstart_one', 'threshold', 'elongation',
        'sstart_mapped_nuc', 'qstart_blastp', 'elongated_length', 'sstart_blastp',
        'qstart_tfastx', 'qstart_prot_diff', 'qstart_nuc_diff',
        'qstart_nuc_elong_right', 'qstart_nuc_short', 'workdir', 'recovered_elongation',
        'raw_recovered_elongation', 'category', 'analysable', "atg_codons_upstream", "atg_explains"
    ]

    schema = {
        'qseqid': pl.Utf8, 'sseqid': pl.Utf8, 'gaps_query': pl.Int32,
        'gaps_subject': pl.Int32, 'qstart_nuc_elong': pl.Int32,
        'sstart_nuc_gt_elongation': pl.Boolean, 'stop_inframe': pl.Boolean,
        'stop_inframe_position': pl.Int32, 'start_inframe': pl.Boolean,
        'start_inframe_position': pl.Int32, 'atg_on_elongated_subject_facing_query_start': pl.Boolean,
        'meth_on_query_facing_subject_start': pl.Boolean, 'sstart_one': pl.Boolean,
        'threshold': pl.Int32, 'elongation': pl.Int32, 'sstart_mapped_nuc': pl.Int32,
        'qstart_blastp': pl.Int32, 'elongated_length': pl.Int32, 'sstart_blastp': pl.Int32,
        'qstart_tfastx': pl.Int32, 'qstart_prot_diff': pl.Int32,
        'qstart_nuc_diff': pl.Int32, 'qstart_nuc_elong_right': pl.Int32, 'qstart_nuc_short': pl.Int32,
        'workdir': pl.Utf8, 'recovered_elongation': pl.Float32, 'raw_recovered_elongation': pl.Int32,
        'category': pl.Utf8, 'analysable': pl.Boolean, "atg_codons_upstream": pl.Int32, "atg_explains": pl.Float32
    }

    (
        results
        .join(thresholds_df, on=["sseqid"], how="inner")
        .join(diff, on=["sseqid"], how="inner")
        .with_columns(
            workdir=pl.lit(workdir),
            recovered_elongation=(((pl.col("qstart_blastp") - 1) * 3 + 1) - pl.col("qstart_nuc_elong")) / ((pl.col("qstart_blastp") - 1) * 3),
            raw_recovered_elongation=((pl.col("qstart_blastp") - 1) * 3 + 1) - pl.col("qstart_nuc_elong"),
            category=pl.lit("short")
        )
        .select(columns)
        .cast(schema)
    ).write_csv(f"{qseqid}_final_res.tsv", separator="\t", include_header=True)

    log.close()
    output.close()


if __name__ == "__main__":
    main()
