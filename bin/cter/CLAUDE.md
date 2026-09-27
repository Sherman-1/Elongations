
# C-ter Extension Analysis

## What this folder does

Scripts in this folder search for C-terminal extension traces in the
non-coding or coding regions of phylogenetic neighbor genomes.

Input: parquet files from parse_blastp_cter.py containing query-subject
pairs where the C-ter of the focal protein is not seen in the homolog.

## Key biological difference with N-ter

The C-ter case is simpler than N-ter because there can only be ONE stop
codon per reading frame. The ribosome stops at the first one. Therefore:

- No need to search for in-frame stop codons
- No need for gap analysis around codons
- No equivalent of the Met-facing-subject-start check
- No equivalent of the ATG-on-elongated-subject-facing-query-start check

We only need to measure how far the nucleotide alignment extends into
the downstream non-coding region of the homolog.

## Big case (slend >= 5)

The homolog has amino acids after the alignment end. The extension trace
might exist in the coding sequence of the homolog but was missed by the
protein alignment.

### big_nuc_search.sh — STATUS: mostly done

1. Converts input parquet to TSV via duckdb
2. Extracts query CDS nucleotide sequence from GFF+FNA
3. For each subject: extracts its CDS nucleotide sequence
4. Runs ssearch36 nucleotide alignment (full query vs all subjects)
5. Outputs: ssearch TSV + qends.tsv

Input parquet columns: qseqid, sseqid, qlend, slend, qend, send
qends.tsv columns: sseqid, qend, qend_mapped (= qend * 3), qlend

Known fixes applied:

- qlend added to qends.tsv (was missing initially)
- is_integer check should include qlend

### parse_big_table.py — STATUS: needs fixes

Computes:

- qend_difference = qend_nuc - qend_mapped
  (positive = nucleotide alignment goes further than protein alignment)
- recovered_elongation = (qend_nuc - qend_mapped) / (qlend * 3)
- raw_recovered_elongation = qend_nuc - qend_mapped

Known issues to fix:

- Formula was incorrect (used wrong column references)
- qend_prot type in schema was pl.Utf8, should be pl.Int32
- .select(columns).cast(schema) was missing before write_csv
- Final .select() was excluding recovered_elongation, raw_recovered_elongation, category

## Small case (slend < 5)

The homolog ends near the end of its coding sequence. The extension trace
is potentially in the DOWNSTREAM non-coding region of the homolog's gene.

### small_nuc_search.sh — STATUS: needs fixes

Same logic as N-ter small but mirrored:

- Elongation direction: downstream (3' side) instead of upstream (5' side)
- Query fragment: last qlend*1.5 nucleotides instead of first qstart*1.5
- Bedtools slop: right side for strand +, left side for strand -

For each subject in the parquet:

1. Compute elongation = qlend_mapped_nuc * 1.5, rounded up to multiple of 3
   - qlend_mapped_nuc = qlend * 3
   - slend = slen - send
   - slend_mapped_nuc = slend * 3 (NOT slen * 3)
2. Elongate the subject gene downstream using bedtools slop (strand-aware):
   - Strand +: slop right on last exon
   - Strand -: slop left on first exon (sorted by genomic position)
3. Extract query C-ter fragment: last qlend_mapped_nuc * 1.5 nucleotides
   - seqkit subseq -r -N:-1 (negative index for C-ter end)
4. Run 3 fragment alignments per subject
5. Record threshold and metrics in subjects_thresholds.tsv

After the loop: run 3 complete alignments (same as N-ter but with full sequences)

Input parquet columns: qseqid, sseqid, qend, send, downstream, Strand, slen, qlen, qlend, slend

Known fixes applied:

- slend_mapped_nuc was incorrectly computed as slen * 3, fixed to (slen - send) * 3
- elongation logged before being calculated — moved after calculation
- is_integer check extended to include qlend and downstream
- threshold column dropped, replaced by elongated_length. threshold only ever
  served N-ter's codon-aligned gap counting (check_atg_gaps), which C-ter does
  not need; it was never read by any C-ter computation. elongated_length is what
  the parser actually needs, and the script computed it without exporting it.

subjects_thresholds.tsv columns:
subject, elongated_length, elongation, slend_mapped_nuc, qend, qlend, slen

The annotated CDS ends at (elongated_length - elongation) in the elongated
subject. This is the anchor the parser needs: an alignment reaching past it has
found the extension trace in the downstream non-coding region.

### parse_small_table.py — STATUS: incomplete

Takes 6 arguments: ssearch_elong.tsv, ssearch_elong.aln, ssearch_short.tsv,
tfastx.tsv, thresholds.tsv, elongated_subjects.faa

#### What is done

- read_thresholds: reads subjects_thresholds.tsv correctly
- parse_ssearch_output: parses alignments, extracts qend for cter side
- Reading alignment tables (table_nuc_elong, table_nuc_short, table_tfastx)
- Diff metrics computation (qend_prot_diff, qend_nuc_diff)
- recovered_elongation formula:
  (qend_nuc_elong - ((qend_blastp - 1) * 3 + 1)) / (qlend_blastp * 3)

#### What is missing

1. No loop over alignments — alignments are parsed but never used.
   Need a loop similar to N-ter that for each subject:

   - Extracts qend_nuc_elong and send_nuc_elong from tables
   - Checks if send_nuc_elong is beyond the elongated region (non-analysable case)
   - Computes recovered_elongation per subject
   - Stores results in a list of dicts
2. No coherence verification between table_nuc_elong and alignments:
   if qend_nuc_elong != alignments[subject]["qend"]: raise ValueError(...)
3. No raise on subject mismatch (currently just logs)
4. qseqid missing from final columns
5. Typo: "recover_elongation" in columns should be "recovered_elongation"
6. .select(columns).cast(schema) missing before write_csv
7. elongated_proteins (SeqIO) read but never used — can be removed
   (no codon checks needed for C-ter)
8. output.txt opened but never used — remove

#### What is NOT needed (unlike N-ter)

- check_atg_gaps function
- check_inframe_codon function
- ATG search in elongated region
- Stop codon search in-frame
- Methionine check in query
- Reading query.faa (no protein sequence analysis needed)
- Reading elongated_subjects.fna nucleotide sequences (no codon analysis)
