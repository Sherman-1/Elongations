# N-ter Extension Analysis

## What this folder does

Scripts in this folder search for N-terminal extension traces in the
non-coding or coding regions of phylogenetic neighbor genomes.

Input: parquet files from parse_blastp_nter.py containing query-subject
pairs where the N-ter of the focal protein is not seen in the homolog.

## Big case (sstart >= 5)

The homolog has amino acids before the alignment start. The extension
trace might exist in the coding sequence of the homolog but was missed
by the protein alignment because the sequences diverged too much.

### big_nuc_search.sh

1. Converts input parquet to TSV via duckdb
2. Extracts query CDS nucleotide sequence from GFF+FNA (gffread)
3. For each subject: extracts its CDS nucleotide sequence
4. Runs ssearch36 nucleotide alignment (full query vs all subjects)
5. Outputs: ssearch TSV + qstarts.tsv

Input parquet columns: qseqid, sseqid, qstart, sstart
qstarts.tsv columns: sseqid, qstart, qstart_mapped (= qstart * 3)

### parse_big_table.py

Computes:

- qstart_difference = qstart_mapped - qstart_nuc
  (positive = nucleotide alignment goes further than protein alignment)
- recovered_elongation = ((qstart-1)*3+1 - qstart_nuc) / ((qstart-1)*3)
- raw_recovered_elongation = (qstart-1)*3+1 - qstart_nuc

## Small case (sstart < 5)

The homolog starts at the very beginning of its coding sequence. The
extension trace is potentially in the UPSTREAM non-coding region of the
homolog's gene. This is the most interesting case biologically.

### small_nuc_search.sh

For each subject in the parquet:

1. Compute elongation = qstart_mapped_nuc * 1.5, rounded up to multiple of 3
   - qstart_mapped_nuc = (qstart - 1) * 3
   - NOTE: upstream value is intentionally NOT used currently (parasitizes results)
2. Elongate the subject gene upstream using bedtools slop (strand-aware):
   - Strand +: slop left on first exon
   - Strand -: slop right on last exon (sorted by genomic position)
   - Middle exons are never touched
3. Extract query N-ter fragment: first qstart_mapped_nuc * 1.5 nucleotides
4. Run 3 fragment alignments per subject (for testing):
   - ssearch nucleotide: fragment vs elongated subject
   - ssearch nucleotide: fragment vs standard subject
   - ssearch protein: fragment translated vs elongated subject translated
5. Record threshold and metrics in subjects_thresholds.tsv

After the loop:
6. Translate full query and all elongated subjects (faTrans -stop)
7. Run 3 complete alignments:

- ssearch nucleotide: full query vs all elongated subjects (+ alignment output)
- ssearch nucleotide: full query vs all standard subjects
- ssearch protein: full query protein vs all elongated subjects protein

Key variables:

- threshold = elongation - 2 → nucleotide just before the annotated ATG in elongated sequence
- elongated_length = actual length of elongated subject (from faSize)

Input parquet columns: qseqid, sseqid, qstart, sstart, upstream, Strand, slen

### parse_small_table.py

Takes 7 arguments: ssearch_elong.tsv, ssearch_elong.aln, ssearch_short.tsv,
tfastx.tsv, thresholds.tsv, elongated_subjects.fna, query.faa

#### Key functions

read_thresholds: reads subjects_thresholds.tsv (headerless)

- Schema: subject, threshold, elongation, sstart_mapped_nuc, qstart, elongated_length, sstart

parse_ssearch_output: parses ssearch -m FA alignment format

- Keeps best HSP per subject (highest bitscore)
- For N-ter: extracts qstart (minimum query position in best HSP)

check_atg_gaps: counts gaps in alignment up to the threshold position

- gaps_query: gaps in the query sequence (shifts codon positions)
- gaps_subject: gaps in the subject sequence

check_inframe_codon: searches for specific codons (ATG or stop) in-frame

- Iterates from start to stop position with step -3 (going upstream)
- IMPORTANT: start must be > stop (assert this)
- Protect against i < 3 to avoid negative indexing

#### Main loop logic

For each subject, two cases:

Case 1 — Non-analysable:
If sstart_nuc_elong > elongation + (sstart_blastp - 1) * 3
→ alignment missed the elongated region entirely
→ set analysable = False, sentinel values for all checks

Case 2 — Normal analysis:

1. check_atg_gaps → count gaps in upstream region
2. Check ATG on elongated subject facing query start:
   - subject_alt_start_position = elongation - ((qstart_blastp-1)*3 - (sstart_blastp-1)*3) - gaps_query + 1
   - Check if codon at that position is ATG
3. Verify the annotated ATG is where expected:
   - putative_start_codon at position (elongation-2+3) must be ATG (assert)
4. Search for stop codon in-frame between threshold and sstart_nuc_elong
   - If found → blocks the elongation hypothesis
5. Search for ATG in-frame in the same region
6. Check Methionine in query facing subject start:
   - Only checked when sstart_blastp == 1 (subject starts at position 1)
   - query[qstart_blastp - 1] == "M"

IMPORTANT: results.append(buffer_dict) must be AFTER buffer_dict.update(),
not before. The dict from check_atg_gaps only has 5 keys, the update adds
the rest.

#### Diff metrics

- qstart_prot_diff = qstart_blastp - qstart_tfastx
  (does protein alignment on elongated subject go further than original blastp?)
- qstart_nuc_diff = qstart_nuc_short - qstart_nuc_elong
  (does nucleotide alignment go further with elongated vs standard subject?)

#### recovered_elongation

recovered_elongation = ((qstart_blastp-1)*3+1 - qstart_nuc_elong) / ((qstart_blastp-1)*3)
raw_recovered_elongation = (qstart_blastp-1)*3+1 - qstart_nuc_elong
