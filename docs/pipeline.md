# Pipeline Data Flow

## Overview

preprocess_input_data → NR (Diamond) → search_homologs → big/small subworkflows

## Step 1 — preprocess_input_data (A1_preprocess.nf)

Inputs:

- input/focal/*.gff + *.fna (focal species)
- input/neighbors/*.gff + *.fna (neighbor species)

Processes:

- concat_gff: concatenate all GFFs into one file
- concat_fasta: concatenate all FNAs into one file
- faidx: index the concatenated FNA (samtools faidx)
- getGeneCoords: extract coding coordinates from GFF (upstream, downstream, strand)
- extract_focal_proteome: extract proteins from focal GFF+FNA using gffread
- extract_neighbors_proteome: extract proteins from all neighbor GFF+FNA
- concat_prot: concatenate all neighbor proteins into one file
- blast_makedb: create a BLAST database from concatenated neighbor proteins

Outputs:

- full_gff: concatenated GFF
- full_fna: concatenated FNA
- index: FNA index (.fai)
- focal_proteome: focal species proteins (.faa)
- local_db: BLAST database of neighbor proteins
- geneCoords: TSV with coding_RNA, upstream, downstream, Strand columns
- local_proteome: concatenated neighbor proteins

## Step 2 — NR (A2_NR.nf)

Inputs:

- focal_proteome
- strain2species dictionary
- eukaryotes dictionary

Processes:

- big_diamond: align focal proteome vs NR database
  - Filters: qcovhsp > 60, scovhsp > 60, ppos > 70
  - Output columns: qseqid, sseqid, qlen, qstart, qlend (= qlen - qend), staxids
- split_ids + split_df: split results into chunks for parallel processing
- parseDiamond (bin/parse_diamond.py): for each protein, count species coverage
  - nter candidate if <= 5% species have qstart <= 15
  - cter candidate if <= 5% species have qlend <= 15

Outputs:

- nter_candidates_IDs: list of protein IDs (no header)
- cter_candidates_IDs: list of protein IDs (no header)

## Step 3 — search_homologs (A3_search_extensions.nf)

### Common to N-ter and C-ter

Inputs:

- nter/cter_candidates_IDs
- local_db, focal_proteome, full_gff, full_fna, index, geneCoords, local_proteome

Processes:

- seqkitGrep: extract candidate protein sequences from focal proteome
- ssearch (or blastp_strict): align candidates against local neighbor proteome
- parseNter / parseCter (bin/parse_blastp_nter.py / bin/parse_blastp_cter.py):
  - Filter for homology (evalue <= 1e-5, qcovhsp >= 70, ppos > 50)
  - Remove bad species (those already covering the extension)
  - Split into big and small based on sstart (N-ter) or slend (C-ter)
  - Output: one parquet file per query, per category (big/small)

### N-ter big subworkflow (nter_big_sstart)

Processes:

- ssearch_big_sstart (bin/nter/big_nuc_search.sh):
  - Input: parquet with columns qseqid, sseqid, qstart, sstart
  - Extracts CDS nucleotide sequences for query and each subject
  - Runs ssearch36 nucleotide alignment (query.fna vs subjects.fna)
  - Output: ssearch TSV + qstarts.tsv (sseqid, qstart, qstart_mapped)
- parse_big_table (bin/nter/parse_big_table.py):
  - Computes recovered_elongation = ((qstart-1)*3+1 - qstart_nuc) / ((qstart-1)*3)
  - Output: {qseqid}_final_res.tsv

### N-ter small subworkflow (nter_small_sstart)

Processes:

- small_elongate_and_align (bin/nter/small_nuc_search.sh):
  - Input: parquet with columns qseqid, sseqid, qstart, sstart, upstream, Strand, slen
  - For each subject:
    - Computes elongation = qstart_mapped_nuc * 1.5 (rounded to multiple of 3)
    - Elongates subject gene upstream with bedtools slop (strand-aware)
    - Extracts query N-ter fragment (first qstart*1.5 nucleotides)
    - Runs 3 alignments: ssearch nuc elongated, ssearch nuc standard, ssearch protein
    - Records threshold (= elongation - 2) and other metrics
  - After loop: runs complete alignments (full query vs all elongated subjects)
  - Output: 11 files (TSV tables, alignments, thresholds, sequences)
- parse_small_table (bin/nter/parse_small_table.py):
  - Parses ssearch alignment output to extract positions and gaps
  - Checks for ATG in-frame in upstream region
  - Checks for stop codon in-frame between putative ATG and annotated start
  - Checks for Met in query facing subject start
  - Computes recovered_elongation
  - Output: {qseqid}_final_res.tsv

### C-ter big subworkflow (cter_big_send)

Same logic as N-ter big but with qend/send instead of qstart/sstart.

- Input parquet columns: qseqid, sseqid, qlend, slend, qend, send
- qends.tsv has columns: sseqid, qend, qend_mapped, qlend
- recovered_elongation = (qend_nuc - qend_mapped) / (qlend * 3)

### C-ter small subworkflow (cter_small_send)

Same elongation logic as N-ter small but downstream instead of upstream.
No ATG/stop codon checks needed (only one stop per reading frame).

- Input parquet columns: qseqid, sseqid, qend, send, downstream, Strand, slen, qlen, qlend, slend
- Elongation direction: downstream (bedtools slop on 3' side, strand-aware)
- Query fragment: last qlend*1.5 nucleotides (seqkit subseq -r -N:-1)
- STATUS: in development, parse_small_table.py incomplete

## Final outputs

All final results are collected with Nextflow collectFile into:

- output/nter/bigFinal
- output/nter/smallFinal
- output/cter/bigFinal (not yet active)
- output/cter/smallFinal (not yet active)
- output/nter/full.blast and output/cter/full.blast (filtered blast results)
