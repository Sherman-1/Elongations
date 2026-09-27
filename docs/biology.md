# Biological Context

## Research question

Do evolutionary events exist where a start or stop codon mutation produces
the elongation of a protein at its N-terminal or C-terminal end?
If so, how can we detect these events computationally?

## Approach

We use comparative genomics between a focal species and its closest phylogenetic
neighbors. The idea is:

1. If a protein in the focal species has a region (N-ter or C-ter) that is
   NOT SEEN in the majority of homologs across eukaryotes, that region is
   evolutionarily recent.
2. To confirm this, we check whether the sequence corresponding to that
   extension exists in neighbor genomes but in a NON-CODING region
   (upstream or downstream of the annotated gene).

This would mean the extension exists in the genome but is not annotated
as coding — evidence of a recent elongation event via codon mutation.

## Pipeline logic step by step

### Step 1 — Diamond vs NR (Non-Redundant database)

We align all proteins of the focal species against the NCBI NR database
using Diamond blastp, filtered on Eukaryota (taxid 2759).

For each protein, we count how many unique eukaryotic species have alignments
that cover the N-terminal (qstart <= 15) or C-terminal (qlend <= 15) region.

If less than 5% of species cover the N-ter → N-ter candidate.
If less than 5% of species cover the C-ter → C-ter candidate.

qlend is computed as qlen - qend in the Diamond output (distance from
alignment end to protein end).

### Step 2 — Local blast against neighbors

Candidates are blasted against the proteomes of phylogenetic neighbors.
Results are filtered for homology (evalue <= 1e-5, qcovhsp >= 70, ppos > 50).

Bad species filtering: for each candidate, species where the extension IS
already seen in at least one protein (qstart < 20 for N-ter, qlend < 20
for C-ter) are removed. These species cannot serve as evidence of elongation
because they already cover the extension.

### Step 3 — Big vs Small split

After filtering, homologs are split based on the position of the alignment
in the SUBJECT (homolog) sequence.

For N-ter, the split is on sstart:

- big (sstart >= 5): the homolog has amino acids before the alignment start.
  The extension trace might be in the CODING sequence of the homolog.
- small (sstart < 5): the homolog starts at the very beginning of its sequence.
  The extension trace is potentially in the UPSTREAM non-coding region.

For C-ter, the split is on slend (= slen - send):

- big (slend >= 5): the homolog has amino acids after the alignment end.
  The extension trace might be in the CODING sequence of the homolog.
- small (slend < 5): the homolog ends near the end of its sequence.
  The extension trace is potentially in the DOWNSTREAM non-coding region.

### Step 4 — Nucleotide search

#### Big case (both N-ter and C-ter)

The protein sequences diverge at their extremities. We ask: is there
nucleotide homology that the protein alignment missed? Nucleotides evolve
slower than amino acids, so we might find a deeper homology signal.

We extract the full CDS nucleotide sequences of query and subject, align
them with ssearch36 nucleotide mode, and check whether the alignment extends
further than the protein alignment did (qstart_nuc < qstart_mapped for N-ter,
qend_nuc > qend_mapped for C-ter).

recovered_elongation measures the fraction of the extension recovered.

#### Small case — N-ter

The homolog starts at the beginning of its coding sequence, so the extension
must be in the upstream non-coding region. We:

1. Elongate the subject gene upstream using bedtools slop (strand-aware).
   Elongation = qstart_mapped_nuc * 1.5, rounded to nearest multiple of 3.
2. Extract the N-ter fragment of the query (first qstart*1.5 nucleotides).
3. Align fragment vs elongated subject (ssearch nucleotide + protein modes).
4. In the elongated subject sequence, check for:
   - An ATG in-frame at the position facing the query start → evidence that
     the homolog has a START codon in its upstream that was missed by annotation.
   - A stop codon in-frame between the putative ATG and the annotated start →
     would block the elongation hypothesis.
   - A Methionine in the query at the position facing the subject start →
     vestige of the ancestral common start codon.

These checks are necessary because in N-ter, there can be multiple ATGs
in the same reading frame. We need to distinguish a true elongation from
a simple annotation error (gene annotated from the wrong ATG).

The threshold (= elongation - 2) points to the nucleotide just before the
annotated ATG of the homolog in the elongated sequence. This is the reference
point for all in-frame codon searches.

#### Small case — C-ter

The C-ter case is conceptually simpler than N-ter. By definition, there can
only be ONE stop codon per reading frame — the ribosome stops at the first
one it encounters. Therefore:

- No need to search for in-frame stop codons (unlike ATG search in N-ter)
- No need for gap analysis around codons
- No equivalent of the Met-facing-subject-start check

We elongate the subject gene downstream using bedtools slop (strand-aware),
align the C-ter fragment of the query against the elongated subject, and
measure how far the alignment extends into the downstream region.

The key metric remains recovered_elongation: how much of the expected
C-ter extension is found in the downstream non-coding region of the homolog.

## Interpretation

A high recovered_elongation means:

- The extension of the focal protein IS found in the neighbor genome
- But it is in a non-coding region of the neighbor
- This suggests the neighbor lost the coding capacity (stop gained or
  start lost) while the focal species retained it

For N-ter specifically, finding an ATG in-frame upstream of the annotated
start in the homolog strengthens the case: the ancestral gene was longer,
the homolog lost its original start codon, and the focal species kept it.
