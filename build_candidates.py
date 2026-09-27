#!/usr/bin/env python3
"""
build_candidates.py - Resolve every genomic question once, and write tables that
downstream work can consume without ever reopening a GFF, a work dir or a
sequence file again.

After this script, the remaining questions are representational, algorithmic or
phylogenetic. Nothing downstream should need to know what a work dir is.

It ANNOTATES, it does not FILTER. Thresholds belong downstream, so they can be
changed without recomputing anything. The two thresholds that do appear here
(--min-ext-frac, --max-gap-frac) define the vocabulary rather than a selection,
they are propagated to every consumer inside this script, and the values used
are written to build_params.tsv next to the tables.

Outputs (in --outdir):

  pairs.tsv         one row per query x subject, every metric, both branches,
                    all three categories (bad / big / small)
  states.tsv        one row per gene x species, a single conservative state;
                    this is what the matrix and the phylogenetic test consume
  genes.tsv         one row per focal LOCUS (isoforms collapsed), for selection
  build_params.tsv  the thresholds and options this run was built with

Columns worth knowing about:

  category           bad = the species already has the extension in an annotated
                     protein; big = the homolog has coding sequence past the
                     alignment; small = the homolog stops at its own end, so any
                     trace can only be non-coding
  ext_frac           fraction of the focal extension covered by the single
                     retained HSP
  spans_boundary     that same HSP covers both the extension and part of the
                     conserved core: the trace is continuous with the gene
                     rather than a fragment floating in the flank
  ext_in_terminal_exon  the focal extension fits inside the first (N-ter) or
                     last (C-ter) CDS exon. When False the extension crosses an
                     exon junction, which points at a splicing event rather than
                     a start/stop mutation
  audit_verdict      independent rescan of the upstream ATG/stop call, done from
                     the sequence rather than from the pipeline's own record.
                     Only ever computed on rows the pipeline itself declared
                     analysable — auditing an unanalysable row compares a real
                     verdict to a missing one and manufactures disagreements
  not_analysable     N-ter row the pipeline could not decide on. Folded into
                     'no_trace' by default; see --separate-not-analysable

Usage:
    python3 build_candidates.py --side nter
    python3 build_candidates.py --side both --outdir analysis/tables
    python3 build_candidates.py --side both --neighbor-gff input/neighbors.gff
"""

import argparse
import json
import re
import sys
from pathlib import Path

import polars as pl

STOPS = {"TAA", "TAG", "TGA"}

# Column names that differ between branches; everything written out is
# normalised so that consumers never need to know which side they are on.
SIDE_COLS = {
    "nter": {"nuc_diff": "qstart_nuc_diff", "ext_len": "qstart_blastp",
             "hsp_lo": "qstart_nuc_elong", "hsp_hi": "qend_nuc_elong"},
    "cter": {"nuc_diff": "qend_nuc_diff", "ext_len": "qlend_blastp",
             "hsp_lo": "qstart_nuc_elong", "hsp_hi": "qend_nuc_elong"},
}

# Columns carried over verbatim from smallFinal when present.
PASSTHROUGH = (
    "elongation", "elongated_length", "analysable", "stop_inframe",
    "start_inframe", "atg_explains", "atg_codons_upstream",
    "recovered_elongation", "raw_recovered_elongation",
    "meth_on_query_facing_subject_start", "sstart_one",
    # conserved ATG facing the focal start: if every neighbour has one, the
    # start is shared and only the annotation differs, which points at the
    # FOCAL annotation being off
    "atg_on_elongated_subject_facing_query_start",
    # kept so downstream scripts can find the work dir without globbing
    "workdir",
)


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--side", choices=["nter", "cter", "both"], default="both")
    p.add_argument("--output", type=Path, default=Path("output"))
    p.add_argument("--work", type=Path, default=Path("work"))
    p.add_argument("--focal-gff", type=Path, default=None,
                   help="focal GFF (default: the single .gff in input/focal/)")
    p.add_argument("--neighbor-gff", type=Path, default=None,
                   help="concatenated neighbour GFF. Only needed to fill "
                        "subject_n_cds_exons; the column is omitted otherwise "
                        "rather than written empty")
    p.add_argument("--outdir", type=Path, default=Path("analysis/tables"))
    p.add_argument("--no-audit", action="store_true",
                   help="skip the independent ATG/stop rescan (faster)")
    p.add_argument("--min-ext-frac", type=float, default=0.8,
                   help="fraction of the extension the single retained HSP must "
                        "cover for the verdict to be 'elongation'")
    p.add_argument("--max-gap-frac", type=float, default=0.2,
                   help="maximum gaps / align_length of that HSP")
    p.add_argument("--coding-flank-priority", choices=["low", "high"],
                   default="low",
                   help="where 'coding_flank' sits when aggregating a species' "
                        "subjects. 'low' (default, historical): any other "
                        "verdict wins, so a species with one unscorable subject "
                        "and one scorable one is scored. 'high': one unscorable "
                        "subject makes the whole species unscorable, which is "
                        "how the downstream matrix already treats the state")
    p.add_argument("--separate-not-analysable", action="store_true",
                   help="emit 'not_analysable' instead of folding undecidable "
                        "N-ter rows into 'no_trace'. Changes the state "
                        "vocabulary; consumers group it with absence either way")
    return p.parse_args()


# ---------------------------------------------------------------- GFF

def read_gff(path: Path):
    """Minimal, fast GFF reader: raw TSV plus regex on the attribute field.

    gff3_parser goes through pandas and needs a lot of memory on concatenated
    files; here only a few fields are ever needed.
    """
    df = pl.read_csv(
        path, separator="\t", has_header=False, comment_prefix="#",
        new_columns=["seqid", "source", "type", "start", "end",
                     "score", "strand", "phase", "attr"],
        infer_schema_length=0, truncate_ragged_lines=True,
    )
    return df.with_columns(pl.col("start").cast(pl.Int64, strict=False),
                           pl.col("end").cast(pl.Int64, strict=False))


def cds_exons(gff):
    """CDS features with their parent mRNA, in transcript order per mRNA."""
    return (
        gff.filter(pl.col("type") == "CDS")
        .with_columns(pl.col("attr").str.extract(r"Parent=([^;]+)").alias("mrna"))
        .filter(pl.col("mrna").is_not_null())
        .select(["mrna", "start", "end", "strand"])
    )


def terminal_exon_length(cds, side):
    """Length of the first (N-ter) or last (C-ter) CDS exon, in transcript order.

    Transcript order follows the strand: on '-' the first exon in transcript
    order is the one with the highest genomic coordinate. Only the length of
    that one exon is needed, because the extension is contained iff it is no
    longer than it.
    """
    plus_first = (side == "nter")
    return (
        cds.with_columns((pl.col("end") - pl.col("start") + 1).alias("len"))
        .with_columns(
            pl.when(pl.col("strand") == "+")
              .then(pl.col("start"))
              .otherwise(-pl.col("end")).alias("_ord"))
        .sort("_ord", descending=not plus_first)
        .group_by("mrna", maintain_order=True)
        .agg(pl.col("len").first().alias("terminal_exon_len"),
             pl.len().alias("n_cds_exons"))
    )


def mrna_to_gene(gff):
    """mRNA ID -> (gene ID, gene name).

    Some annotations declare only CDS features (three of the yeast GFFs do).
    An inner join on an empty mRNA table returns nothing without erroring and
    the locus level silently disappears, so fall back to the CDS Parent and say
    so rather than emit a table full of nulls.
    """
    m = (gff.filter(pl.col("type") == "mRNA")
         .with_columns(
             pl.col("attr").str.extract(r"ID=([^;]+)").alias("mrna"),
             pl.col("attr").str.extract(r"Parent=([^;]+)").alias("gene_id"),
             pl.col("attr").str.extract(r"gene=([^;]+)").alias("gene_name"))
         .select(["mrna", "gene_id", "gene_name"])
         .filter(pl.col("mrna").is_not_null()))

    if m.is_empty():
        print("WARNING: no mRNA feature in the focal GFF; deriving the locus "
              "level from the CDS Parent. Isoforms of one gene will NOT be "
              "collapsed.", file=sys.stderr)
        m = (gff.filter(pl.col("type") == "CDS")
             .with_columns(
                 pl.col("attr").str.extract(r"Parent=([^;]+)").alias("mrna"))
             .filter(pl.col("mrna").is_not_null())
             .select("mrna")
             .unique()
             .with_columns(pl.col("mrna").alias("gene_id"),
                           pl.lit(None, dtype=pl.Utf8).alias("gene_name")))

    m = m.unique(subset=["mrna"], keep="first")
    # a locus with no symbol is still a locus; keep the row groupable
    return m.with_columns(
        pl.col("gene_name").fill_null(pl.col("gene_id")).alias("gene_name"))


def count_subject_exons(path: Path):
    """CDS exons per transcript in the concatenated neighbour GFF.

    Streamed line by line: the concatenated file runs to a couple of gigabytes
    and only one counter per Parent is ever needed.
    """
    parent = re.compile(r"(?:^|;)Parent=([^;\t\n]+)")
    counts = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 9)
            if len(f) < 9 or f[2] != "CDS":
                continue
            m = parent.search(f[8])
            if m:
                for pid in m.group(1).split(","):
                    counts[pid] = counts.get(pid, 0) + 1
    return pl.DataFrame({"sseqid": list(counts.keys()),
                         "subject_n_cds_exons": list(counts.values())},
                        schema={"sseqid": pl.Utf8,
                                "subject_n_cds_exons": pl.Int64})


# ---------------------------------------------------------------- work dirs

def locate_work(gene, recorded, work: Path, cache):
    """Directory holding the ssearch tables and sequences for a gene."""
    if gene in cache:
        return cache[gene]
    d = work / recorded if recorded else None
    if d is not None and (d / "elongated_subjects.fna").exists():
        out = (d / "elongated_subjects.fna").resolve().parent
    else:
        # the parse work dir may have been deleted to force a rerun; the
        # sequences survive in the elongate dir, found by its thresholds file
        hits = sorted(work.glob(f"*/*/{gene}_subjects_thresholds.tsv"))
        if len(hits) > 1:
            print(f"  note: {len(hits)} work dirs for {gene}, taking the last "
                  f"({hits[-1].parent})", file=sys.stderr)
        out = hits[-1].parent if hits else None
    cache[gene] = out
    return out


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


def audit_upstream(seq, elongation):
    """Independent rescan: walk in frame from the annotated ATG towards 5'.

    Returns (verdict, codons_to_first_atg, codons_to_first_stop). verdict is
    'annotation_doubt' when a usable ATG is reached before any stop.
    """
    first_atg = first_stop = None
    i, k = elongation - 3, 1
    while i >= 0:
        codon = seq[i:i + 3].upper()
        if codon in STOPS and first_stop is None:
            first_stop = k
        if codon == "ATG" and first_atg is None:
            first_atg = k
        if first_stop and first_atg:
            break
        i -= 3
        k += 1
    usable = first_atg is not None and (first_stop is None or first_atg < first_stop)
    return ("annotation_doubt" if usable else "elongation"), first_atg, first_stop


def is_true(v):
    return str(v).strip().lower() in {"true", "1", "yes", "t"}


# ---------------------------------------------------------------- per side

def load_ssearch(tbl: Path):
    """Best HSP per subject plus the HSP count, as plain dicts.

    Filtering a DataFrame once per subject is quadratic and shows up badly on
    the mammal datasets, where a single gene can carry hundreds of subjects.
    """
    a = pl.read_csv(tbl, separator="\t", infer_schema_length=10000,
                    schema_overrides={"evalue": pl.Float32})
    best = (a.sort("bitscore", descending=True)
             .unique(subset=["sseqid"], keep="first", maintain_order=True))
    n_hsp = {}
    for sid in a["sseqid"].to_list():
        n_hsp[sid] = n_hsp.get(sid, 0) + 1
    by_sid = {r["sseqid"]: r for r in best.iter_rows(named=True)}
    return by_sid, n_hsp


def build_side(side, args, gene_map, term_exon, subject_exons):
    src = args.output / side
    cols = SIDE_COLS[side]

    small = pl.read_csv(src / "smallFinal", separator="\t",
                        infer_schema_length=10000)
    blast = pl.read_csv(src / "full.blast", separator="\t",
                        infer_schema_length=10000)
    big = (pl.read_csv(src / "bigFinal", separator="\t", infer_schema_length=10000)
           if (src / "bigFinal").exists() else None)
    bad = (pl.read_csv(src / "bad_species.tsv", separator="\t")
           if (src / "bad_species.tsv").exists() else None)

    keep = [c for c in ["qseqid", "sseqid", "species", "qlen", "qstart", "qend",
                        "sstart", "send", "slen", "evalue", "bitscore",
                        "qcovhsp", "ppos", "scovhsp"] if c in blast.columns]
    bl = blast.select(keep)
    # one blastp row per pair, chosen explicitly rather than by file order
    if "bitscore" in bl.columns:
        bl = bl.sort("bitscore", descending=True, nulls_last=True)
    elif "evalue" in bl.columns:
        bl = bl.sort("evalue", descending=False, nulls_last=True)
    bl = bl.unique(subset=["qseqid", "sseqid"], keep="first", maintain_order=True)

    rows = []
    cache = {}
    aln_cache = {}
    audited = skipped_audit = 0

    # ---- small: the informative category, needs the work dirs
    for gene in small["qseqid"].unique().to_list():
        sub = small.filter(pl.col("qseqid") == gene)
        wd = locate_work(gene, sub["workdir"][0] if "workdir" in sub.columns
                         else None, args.work, cache)

        best = nhsp = None
        seqs = {}
        if wd is not None:
            tbl = wd / f"{gene}_complete_small_ssearch_elong.tsv"
            if tbl.exists():
                if tbl not in aln_cache:
                    aln_cache[tbl] = load_ssearch(tbl)
                best, nhsp = aln_cache[tbl]
            fa = wd / "elongated_subjects.fna"
            if fa.exists() and not args.no_audit:
                seqs = read_fasta(fa)

        for r in sub.iter_rows(named=True):
            sid = r["sseqid"]
            rec = {"side": side, "qseqid": gene, "sseqid": sid,
                   "category": "small"}

            ext_aa = r.get(cols["ext_len"])
            ext_aa = int(ext_aa) if ext_aa is not None else None
            ext_nt = ((ext_aa - 1) * 3 if side == "nter" else ext_aa * 3) \
                if ext_aa is not None else None
            rec["ext_aa"] = (ext_aa - 1) if (side == "nter" and ext_aa) else ext_aa
            rec["ext_nt"] = ext_nt

            for k in PASSTHROUGH:
                if k in small.columns:
                    rec[k] = r.get(k)
            rec["nuc_diff"] = r.get(cols["nuc_diff"])

            # ---- alignment shape, from the raw ssearch table
            b = best.get(sid) if best else None
            if b is not None:
                lo, hi = int(b[cols["hsp_lo"]]), int(b[cols["hsp_hi"]])
                qlen_nuc = int(b["qlen_nuc_elong"])
                if ext_nt:
                    if side == "nter":
                        cov = max(0, min(hi, ext_nt) - lo + 1)
                        spans = lo <= ext_nt < hi
                    else:
                        estart = qlen_nuc - ext_nt + 1
                        cov = max(0, hi - max(lo, estart) + 1)
                        spans = lo < estart <= hi
                else:
                    cov, spans = None, None
                rec.update({
                    "n_hsp": nhsp.get(sid, 1),
                    "align_length": int(b["align_length"]),
                    "gaps": int(b["gaps"]),
                    "pident": float(b["pident"]),
                    "ext_covered": cov,
                    "ext_frac": round(cov / ext_nt, 4)
                                if (cov is not None and ext_nt) else None,
                    "spans_boundary": spans,
                })

            # ---- independent audit of the upstream call (N-ter only)
            #
            # The analysable guard is not an optimisation. Auditing a row the
            # pipeline could not decide on compares a real verdict to a missing
            # one; without it the audit reports ~880 disagreements that are
            # nothing but that mismatch.
            if side == "nter" and seqs and sid in seqs and r.get("elongation"):
                if is_true(r.get("analysable")):
                    v, atg_k, stop_k = audit_upstream(seqs[sid],
                                                      int(r["elongation"]))
                    rec["audit_verdict"] = v
                    rec["audit_atg_codons"] = atg_k
                    rec["audit_stop_codons"] = stop_k
                    audited += 1
                else:
                    skipped_audit += 1

            rows.append(rec)

    if side == "nter" and not args.no_audit:
        print(f"  audit: {audited} rows rescanned, {skipped_audit} skipped as "
              f"not analysable", file=sys.stderr)

    # ---- big: no work dir needed, it is the non-informative category
    if big is not None:
        for r in big.iter_rows(named=True):
            rows.append({"side": side, "qseqid": r["qseqid"],
                         "sseqid": r["sseqid"], "category": "big",
                         "recovered_elongation": r.get("recovered_elongation"),
                         "raw_recovered_elongation":
                             r.get("raw_recovered_elongation")})

    pairs = pl.DataFrame(rows, infer_schema_length=None)

    # ---- attach blastp context and derive the internal-indel proxy
    pairs = pairs.join(bl, on=["qseqid", "sseqid"], how="left")
    if {"qstart", "qend", "sstart", "send"} <= set(pairs.columns):
        pairs = pairs.with_columns(
            ((pl.col("qend") - pl.col("qstart")) -
             (pl.col("send") - pl.col("sstart"))).abs().alias("blastp_indel"))

    # ---- focal locus and exon structure
    pairs = (pairs.join(gene_map, left_on="qseqid", right_on="mrna", how="left")
                  .join(term_exon, left_on="qseqid", right_on="mrna", how="left"))
    if {"terminal_exon_len", "ext_nt"} <= set(pairs.columns):
        pairs = pairs.with_columns(
            (pl.col("ext_nt") <= pl.col("terminal_exon_len"))
            .alias("ext_in_terminal_exon"))

    # subject exon counts come from the NEIGHBOUR annotation. Joining the focal
    # transcript ids against subject ids can only ever match nothing, so the
    # column is filled when the neighbour GFF is supplied and omitted otherwise.
    if subject_exons is not None:
        pairs = pairs.join(subject_exons, on="sseqid", how="left")

    # ---- bad species: one row per (gene, species), no subject
    if bad is not None:
        bad_rows = (bad.with_columns(
                        pl.col("bad_species").fill_null("")
                          .str.split(";"))
                    .explode("bad_species")
                    .rename({"bad_species": "species"})
                    .filter(pl.col("species") != "")
                    .select(["qseqid", "species"])
                    .unique()
                    .with_columns(pl.lit(side).alias("side"),
                                  pl.lit("bad").alias("category")))
        bad_rows = (bad_rows
                    .join(gene_map, left_on="qseqid", right_on="mrna", how="left"))
        pairs = pl.concat([pairs, bad_rows], how="diagonal")

    return pairs


# ---------------------------------------------------------------- states

def verdict_expr(min_ext_frac, max_gap_frac, separate_not_analysable):
    """Per-row verdict.

    Two things are being separated here.

    N-ter can tell mis-annotation from elongation; C-ter cannot, because a stop
    is unique per reading frame while several ATGs can be.

    And a trace is not a trace just because the alignment gained ground. A
    handful of nucleotides recovered, or a stretch broken by a large gap, says
    almost nothing: 'elongation' is reserved for a SINGLE local alignment that
    covers most of the extension AND reaches into the conserved core. Anything
    that gains ground without meeting that bar becomes 'weak_trace' rather than
    being silently promoted or discarded — the distinction stays visible on the
    figures and the thresholds stay adjustable.
    """
    b = lambda c: pl.col(c).cast(pl.Boolean, strict=False).fill_null(False)
    gap_frac = (pl.col("gaps").cast(pl.Float64, strict=False)
                / pl.col("align_length").cast(pl.Float64, strict=False))
    solid = (b("spans_boundary")
             & (pl.col("ext_frac") >= min_ext_frac)
             & (gap_frac.fill_null(0.0) <= max_gap_frac))
    undecided = "not_analysable" if separate_not_analysable else "no_trace"
    return (
        pl.when(pl.col("category") == "bad").then(pl.lit("complete"))
        .when(pl.col("category") == "big").then(pl.lit("coding_flank"))
        .when((pl.col("side") == "nter") & ~b("analysable"))
            .then(pl.lit(undecided))
        .when((pl.col("side") == "cter") & (pl.col("nuc_diff") <= 0))
            .then(pl.lit("no_trace"))
        .when((pl.col("side") == "nter") & b("start_inframe") & ~b("stop_inframe"))
            .then(pl.lit("annotation_doubt"))
        .when(solid).then(pl.lit("elongation"))
        .otherwise(pl.lit("weak_trace"))
    )


def priority_order(coding_flank_priority):
    """Most conservative first: one doubtful subject outweighs any number of
    supporting ones, exactly as one covering transcript condemns a whole species.

    Where 'coding_flank' belongs is a scientific choice, not a technical one.
    Historically it sat last, so a species with one unscorable subject and one
    scorable one was scored on the scorable one — while the downstream matrix
    treats the very same state as disqualifying. Both readings are defensible;
    they must not be held at the same time silently.
    """
    if coding_flank_priority == "high":
        base = ["complete", "coding_flank", "annotation_doubt", "weak_trace",
                "elongation", "not_analysable", "no_trace"]
    else:
        base = ["complete", "annotation_doubt", "weak_trace", "elongation",
                "not_analysable", "no_trace", "coding_flank"]
    return base


def build_states(pairs, priority):
    rank = {v: i for i, v in enumerate(priority)}
    dropped = pairs.filter(pl.col("species").is_null()).height
    if dropped:
        print(f"  {dropped} pairs with no species, dropped from states.tsv",
              file=sys.stderr)
    unknown = set(pairs["verdict"].drop_nulls().to_list()) - set(rank)
    if unknown:
        sys.exit(f"verdicts with no priority defined: {sorted(unknown)}")
    return (
        pairs.filter(pl.col("species").is_not_null())
        .with_columns(pl.col("verdict").replace_strict(rank, default=None)
                        .alias("_r"))
        .group_by(["side", "qseqid", "species"])
        .agg(pl.col("_r").min())
        .with_columns(pl.col("_r").replace_strict({v: k for k, v in rank.items()},
                                                  default="absent")
                        .alias("state"))
        .drop("_r")
        .sort(["side", "qseqid", "species"])
    )


def build_genes(pairs, min_ext_frac):
    """One row per focal LOCUS. Isoforms sharing an extension length are the
    same observation, so the key is (gene_id, ext_aa).

    n_species_contiguous applies the same --min-ext-frac as the verdict does.
    A hard-coded 0.8 here would let the two disagree, and this is the column
    genes.tsv is sorted on.
    """
    sm = pairs.filter(pl.col("category") == "small")
    if sm.is_empty():
        return pl.DataFrame()
    return (
        sm.group_by(["side", "gene_id", "gene_name", "ext_aa"])
        .agg(
            pl.col("qseqid").n_unique().alias("n_isoforms"),
            pl.col("species").n_unique().alias("n_species"),
            pl.col("species").filter(pl.col("verdict") == "elongation")
              .n_unique().alias("n_species_elongation"),
            pl.col("species").filter(
                pl.col("verdict") == "annotation_doubt")
              .n_unique().alias("n_species_doubt"),
            pl.col("species").filter(
                pl.col("spans_boundary").cast(pl.Boolean, strict=False)
                & (pl.col("ext_frac") >= min_ext_frac))
              .n_unique().alias("n_species_contiguous"),
            pl.col("recovered_elongation").median().alias("rec_median"),
            pl.col("ext_in_terminal_exon").cast(pl.Boolean, strict=False)
              .all().alias("ext_in_terminal_exon"),
        )
        .sort(["n_species_contiguous", "n_species_elongation"],
              descending=[True, True])
    )


# ---------------------------------------------------------------- main

def main():
    args = parse_args()
    sides = ["nter", "cter"] if args.side == "both" else [args.side]

    focal = args.focal_gff
    if focal is None:
        hits = sorted(Path("input/focal").glob("*.gff"))
        if not hits:
            sys.exit("no focal GFF found; pass --focal-gff")
        focal = hits[0]
    print(f"focal GFF: {focal}", file=sys.stderr)

    gff = read_gff(focal)
    gene_map = mrna_to_gene(gff)
    cds = cds_exons(gff)

    subject_exons = None
    if args.neighbor_gff:
        if not args.neighbor_gff.exists():
            sys.exit(f"missing {args.neighbor_gff}")
        print(f"neighbour GFF: {args.neighbor_gff} (counting CDS exons)",
              file=sys.stderr)
        subject_exons = count_subject_exons(args.neighbor_gff)
        print(f"  {subject_exons.height} subject transcripts", file=sys.stderr)
    else:
        print("no --neighbor-gff: subject_n_cds_exons will not be written",
              file=sys.stderr)

    args.outdir.mkdir(parents=True, exist_ok=True)
    all_pairs = []

    for side in sides:
        if not (args.output / side / "smallFinal").exists():
            print(f"skipping {side}: no smallFinal", file=sys.stderr)
            continue
        term = terminal_exon_length(cds, side)
        p = build_side(side, args, gene_map,
                       term.select(["mrna", "terminal_exon_len"]), subject_exons)
        all_pairs.append(p)
        print(f"{side}: {p.height} pairs", file=sys.stderr)

    if not all_pairs:
        sys.exit("nothing to build")

    pairs = pl.concat(all_pairs, how="diagonal")

    # verdict_expr references N-ter-only columns; with --side cter alone they
    # never appear, and polars needs them to exist even on branches that can
    # never be true for those rows
    for c, dt in (("analysable", pl.Boolean), ("start_inframe", pl.Boolean),
                  ("stop_inframe", pl.Boolean), ("spans_boundary", pl.Boolean),
                  ("nuc_diff", pl.Int64), ("ext_frac", pl.Float64),
                  ("gaps", pl.Int64), ("align_length", pl.Int64)):
        if c not in pairs.columns:
            pairs = pairs.with_columns(pl.lit(None, dtype=dt).alias(c))

    pairs = pairs.with_columns(
        verdict_expr(args.min_ext_frac, args.max_gap_frac,
                     args.separate_not_analysable).alias("verdict"))
    # traceable even when folded into no_trace
    pairs = pairs.with_columns(
        ((pl.col("side") == "nter")
         & ~pl.col("analysable").cast(pl.Boolean, strict=False).fill_null(False)
         & (pl.col("category") == "small")).alias("not_analysable"))

    priority = priority_order(args.coding_flank_priority)
    states = build_states(pairs, priority)
    genes = build_genes(pairs, args.min_ext_frac)

    pairs.write_csv(args.outdir / "pairs.tsv", separator="\t")
    states.write_csv(args.outdir / "states.tsv", separator="\t")
    if not genes.is_empty():
        genes.write_csv(args.outdir / "genes.tsv", separator="\t")

    params = {
        "min_ext_frac": args.min_ext_frac,
        "max_gap_frac": args.max_gap_frac,
        "coding_flank_priority": args.coding_flank_priority,
        "separate_not_analysable": args.separate_not_analysable,
        "audit": not args.no_audit,
        "focal_gff": str(focal),
        "neighbor_gff": str(args.neighbor_gff) if args.neighbor_gff else "",
        "sides": ",".join(sides),
        "priority": ">".join(priority),
    }
    pl.DataFrame({"param": list(params), "value": [str(v) for v in params.values()]}
                 ).write_csv(args.outdir / "build_params.tsv", separator="\t")

    print()
    print(pairs.group_by(["side", "verdict"]).len().sort(["side", "verdict"]))
    print()
    print(f"pairs : {pairs.height:6d} -> {args.outdir/'pairs.tsv'}")
    print(f"states: {states.height:6d} -> {args.outdir/'states.tsv'}")
    if not genes.is_empty():
        print(f"genes : {genes.height:6d} -> {args.outdir/'genes.tsv'}")
    print(f"params: {args.outdir/'build_params.tsv'}")


if __name__ == "__main__":
    main()