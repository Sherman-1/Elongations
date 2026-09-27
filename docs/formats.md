# File Formats Reference

## Diamond output (after awk filtering in big_diamond)

Columns: qseqid, sseqid, qlen, qstart, qlend, staxids
- qlend = qlen - qend (computed in awk)
- staxids can contain multiple taxids separated by ";"
- Filters applied: qcovhsp > 60, scovhsp > 60, ppos > 70

## parseDiamond output

nter_candidates_*.tsv and cter_candidates_*.tsv:
- Single column, no header
- Contains qseqid values only

statistics_*.tsv:
- Full statistics per protein with header
- Columns include: qseqid, nb_unique_species, nter_q_15, cter_q_15,
  qlen, qstart_median, qstart_std, cter_median, cter_std,
  percent_nter_q_15, percent_cter_q_15

## ssearch / blastp_strict output

Standard blast tabular format with header.
Columns: qseqid, sseqid, qlen, slen, qstart, qend, sstart, send,
evalue, qcovhsp, scovhsp, ppos

## geneCoords (from getGeneCoords)

TSV with header.
Columns: coding_RNA, upstream, downstream, Strand
- upstream: nucleotide distance between gene start and CDS start
- downstream: nucleotide distance between CDS end and gene end
- Strand: "+" or "-"

## parse_blastp_nter.py outputs

### Parquet files (one per query)

{qseqid}_big.parquet:
- Columns: qseqid, sseqid, qstart, sstart

{qseqid}_small.parquet:
- Columns: qseqid, sseqid, qstart, sstart, upstream, Strand, slen

### Other outputs
- full.blast: all filtered blast results (TSV with header)
- filtered.blast: after bad species removal
- bad_species.tsv: qseqid + list of bad species (semicolon-separated)
- mRNA_species.tsv: ID, species (header)

## parse_blastp_cter.py outputs

### Parquet files (one per query)

{qseqid}_big.parquet:
- Columns: qseqid, sseqid, qlend, slend, qend, send

{qseqid}_small.parquet:
- Columns: qseqid, sseqid, qend, send, downstream, Strand, slen, qlen, qlend, slend

### Other outputs
- Same as N-ter (full.blast, filtered.blast, bad_species.tsv, mRNA_species.tsv)

## subjects_thresholds.tsv (produced by small_nuc_search.sh)

Headerless TSV, one row per subject.

### N-ter version
Columns: subject, threshold, elongation, sstart_mapped_nuc, qstart, elongated_length, sstart
- threshold = elongation - 2 (position just before annotated ATG in elongated sequence)
- elongation = qstart_mapped_nuc * 1.5 rounded to multiple of 3
- sstart_mapped_nuc = (sstart - 1) * 3
- elongated_length = length of elongated subject sequence (from faSize)

### C-ter version
Columns: subject, elongated_length, elongation, slend_mapped_nuc, qend, qlend, slen
- elongated_length = length of elongated subject sequence (from faSize)
- elongation = qlend_mapped_nuc * 1.5 rounded to multiple of 3
- slend_mapped_nuc = (slen - send) * 3
- the annotated CDS ends at (elongated_length - elongation); everything past it
  is the downstream non-coding region
- no threshold column: unlike N-ter, C-ter needs no codon-aligned anchor, since
  there is one stop per reading frame and thus no in-frame codon search

## qstarts.tsv / qends.tsv (produced by big_nuc_search.sh)

TSV with header.

### N-ter (qstarts.tsv)
Columns: sseqid, qstart, qstart_mapped
- qstart_mapped = qstart * 3

### C-ter (qends.tsv)
Columns: sseqid, qend, qend_mapped, qlend
- qend_mapped = qend * 3

## ssearch output files (produced by small_nuc_search.sh)

All have the same 15-column format with header:
qseqid, qlen_nuc_*, sseqid, slen_nuc_*, pident, align_length, bar, gaps,
qstart_nuc_*, qend_nuc_*, sstart_nuc_*, send_nuc_*, evalue, bitscore, CIGAR

The suffix indicates the type:
- _elong: query vs elongated subject (nucleotide)
- _short: query vs standard subject (nucleotide)
- _tfastx: query vs elongated subject (protein)

And the prefix indicates scope:
- fragment_: N-ter/C-ter fragment only (per subject, inside loop)
- complete_: full query vs all subjects (after loop)

## final_res.tsv (produced by parse_big_table.py and parse_small_table.py)

### N-ter big
Columns: qseqid, sseqid, qstart_mapped, qstart_nuc, qstart_difference,
qstart, category, recovered_elongation, raw_recovered_elongation, qstart_prot

### N-ter small
Columns: qseqid, sseqid, gaps_query, gaps_subject, qstart_nuc_elong,
sstart_nuc_gt_elongation, stop_inframe, stop_inframe_position,
start_inframe, start_inframe_position, atg_on_elongated_subject_facing_query_start,
meth_on_query_facing_subject_start, sstart_one, threshold, elongation,
sstart_mapped_nuc, qstart_blastp, elongated_length, sstart_blastp,
qstart_tfastx, qstart_prot_diff, qstart_nuc_diff,
qstart_nuc_elong_right, qstart_nuc_short, workdir, recovered_elongation,
raw_recovered_elongation, category, analysable

### C-ter big
Columns: qseqid, sseqid, qend_mapped, qend_nuc, qend_difference,
qend, recovered_elongation, raw_recovered_elongation, category, qend_prot

### C-ter small
STATUS: incomplete — columns not finalized yet