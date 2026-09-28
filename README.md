# Detection and validation of protein terminal extensions by comparative genomics

Nextflow pipeline that identifies, in a focal eukaryotic species, N-terminal and C-terminal protein segments absent from the annotated homologs of other species, and searches for these segments in the unannotated genomic DNA of a set of neighboring species. The pipeline is generic: it applies to any eukaryotic focal species and any set of neighboring species for which a genome assembly and a GFF3 annotation are available.

*Developed by Simon Herman — BIM team.*

## Contents

1. [Background and objective](#1-background-and-objective)
2. [Materials](#2-materials)
3. [Methods](#3-methods)
4. [Outputs](#4-outputs)
5. [Usage](#5-usage)
6. [Repository layout](#6-repository-layout)
7. [References](#7-references)

---

## 1. Background and objective

Automatic genome annotation frequently yields erroneous start or stop codons, causing conserved coding segments to be overlooked. When a well-annotated species (the *focal* species) is compared with less well-annotated neighbors, a number of its proteins appear longer than their annotated homologs. For each such length difference, three hypotheses can be considered:

1. a **genuine innovation** of the focal lineage;
2. an **annotation error in the focal species** (start codon placed too far upstream, or stop codon too far downstream);
3. an **under-annotation of the neighboring species**, the segment being present in its genome but not predicted as coding.

The pipeline aims to discriminate between these hypotheses in three stages: (i) identification, by searching the NR database, of focal proteins whose terminus is found in virtually no homolog; (ii) alignment of these candidates against the proteomes of the neighboring species; (iii) for each candidate/homolog pair, virtual extension of the neighbor's gene into genomic DNA and nucleotide alignment to test for the presence of the missing segment.

---

## 2. Materials

### 2.1 Required input data and constraints

The pipeline takes as input:

- a eukaryotic **focal species**: genome assembly (`.fna`) and GFF3 annotation (`.gff`) in `params.focal_dir`;
- a set of **neighboring species**, of any size: for each, a genome (`[Species].fna`) and a GFF3 annotation (`[Species].gff`) in `params.neighbors_dir`.

A single constraint applies. Because the genomes and annotations of all species are concatenated into one indexed sequence set (§3.2), **sequence identifiers (chromosomes, scaffolds, contigs) must be unique across species**: identical between the `.fna` and `.gff` of a given species, but distinct from those of every other species. The `checkAnnot` module verifies this uniqueness before concatenation.

An optional Newick tree describing the phylogenetic relationships among the species may be provided for downstream interpretation of the results.

### 2.2 External databases

- **NCBI NR** in DIAMOND format (`nr_2.0.13.dmnd`; database format version 3, DIAMOND build 151; 707,028,945 sequences, 2.73 × 10¹¹ residues), queried with a restriction to Eukaryota (taxid 2759).
- **NCBI Taxonomy** (`taxdump`), used for strain → species mapping and for defining the set of eukaryotic taxids.

### 2.3 Software

All executable paths are configurable in `nextflow.config` (see §5.2). Versions below are those used during development, as recorded by `get_versions.sh` (full output in `versions.md`).

| Tool | Version | Role in the pipeline | Reference |
| --- | --- | --- | --- |
| Nextflow (DSL2) | 25.10.4 (Java 21.0.11) | Workflow orchestration | Di Tommaso et al., 2017 |
| Singularity | 1.4.5 | Containerization (AGAT) | Kurtzer et al., 2017 |
| PBS Pro | *to be specified* | HPC scheduling | — |
| DIAMOND | 2.1.8 | Similarity search against NR | Buchfink et al., 2021 |
| BLAST+ (`makeblastdb`, `blastp`) | 2.16.0 | Local protein database, protein alignments | Camacho et al., 2009 |
| FASTA36 (`ssearch36`, `tfasty36`) | 36.3.8i | Smith-Waterman alignments: protein/protein, nucleotide/nucleotide, protein/translated DNA | Pearson & Lipman, 1988; Pearson, 1991 |
| gffread | 0.12.7 | Extraction of CDS and protein sequences from GFF + genome | Pertea & Pertea, 2020 |
| BEDTools (`slop`) | 2.31.1 | Extension of genomic coordinates | Quinlan & Hall, 2010 |
| SAMtools (`faidx`) | 1.21 | Genome indexing | Danecek et al., 2021 |
| SeqKit | 2.6.1 | FASTA manipulation | Shen et al., 2016 |
| TaxonKit | 0.20.0 | NCBI taxonomy handling | Shen & Ren, 2021 |
| AGAT | *to be specified* (`agat.sif` container) | GFF standardization | Dainat |
| UCSC utilities (`faTrans`, `faSize`) | binaries dated 2024-02-09 (no version flag) | DNA → protein translation, sequence sizes | Kent et al., 2002 |
| DuckDB | 1.2.0 | Parquet → TSV conversion | Raasveldt & Mühleisen, 2019 |
| Python | 3.13.5 | Table parsing and processing | — |
| Polars | 1.37.1 | Dataframe processing | — |
| Biopython | 1.85 | Sequence manipulation | Cock et al., 2009 |
| gff3_parser | 0.0.5 | GFF3 parsing | — |
| g++ | 14.2.0 | Compilation of `src/` | — |

---

## 3. Methods

### 3.1 Overview

The pipeline consists of four stages (Figure 1), identical regardless of the species analyzed: data preprocessing (§3.2), candidate identification by searching NR (§3.3), homolog search in the neighboring species (§3.4), and nucleotide-level analysis of the extensions (§3.5). The N-terminal ("nter") and C-terminal ("cter") ends are handled by two symmetric branches.

```
                              main.nf → align
                                    │
             ┌──────────────────────┼──────────────────────┐
             ▼                      ▼                      ▼
     A1_preprocess          create_taxon_maps            A2_NR
  proteome extraction       strain2species        DIAMOND vs NR (euk.)
  concat / faidx / makedb   eukaryote taxids      parseDiamond → candidates
  getGeneCoords                                   nter / cter
             │                                              │
             └──────────────────────┬───────────────────────┘
                                    ▼
                         A3_search_extensions
                  ssearch36 candidates vs local proteome
                  parseNter / parseCter → big | small
                                    │
                     ┌──────────────┴──────────────┐
                     ▼                             ▼
             big_sstart / big_send       small_sstart / small_send
             nucleotide ssearch (CDS)    elongation (bedtools slop)
             + parse                     + multiple alignments
                                         + parse (gaps, codons, recovered elongation)
                                    │
                                    ▼
                          output/nter/  output/cter/
```

*Figure 1. Pipeline architecture.*

### 3.2 Data preprocessing (`A1_preprocess.nf`, `create_taxon_maps.nf`)

**Proteome extraction.** For each species, the proteome is derived from the annotation and genome with `gffread` (option `-J`, which discards transcripts lacking a start codon or containing an in-frame stop), and gaps are removed with `seqkit`:

```bash
gffread -J -y - -g ${fna} ${gff} | seqkit seq --remove-gaps - > ${species}.faa
```

**Concatenated datasets.** The GFF files, genomes and proteomes of all neighboring species are concatenated (`full_gff`, `full_fna`, `local_proteome`). The concatenated genome is indexed with `samtools faidx`, and a local protein database is built with `makeblastdb` (`-dbtype prot -parse_seqids -hash_index`).

**Gene coordinates.** For each gene, the `getGeneCoords` module computes, in a strand-aware manner, the distance between the CDS start and the gene start (*upstream* region) and between the CDS end and the gene end (*downstream* region). These distances bound the elongation applied in §3.5.1.

**Taxonomy tables.** From the NCBI `taxdump`, `taxonkit` produces a strain → species mapping table (`strain2species.csv`) and the list of eukaryotic taxids descending from taxid 2759 (`eukaryotes.csv`).

### 3.3 Candidate identification by searching NR (`A2_NR.nf`)

**Similarity search.** The focal proteome is aligned against NR with `diamond blastp` in `--fast` mode, restricted to Eukaryota, with no limit on the number of targets:

```bash
diamond blastp --query ${query} --db ${db} --taxonlist 2759 \
  --outfmt 6 qseqid sseqid qlen qstart qend qcovhsp scovhsp ppos staxids \
  --max-target-seqs 0 --fast -e 0.00001
```

Alignments are then filtered on query coverage (`qcovhsp` > 60%), subject coverage (`scovhsp` > 60%) and percentage of positive positions (`ppos` > 70%).

**Parallelization.** Query identifiers are distributed into `n` groups (30 by default; `split_ids`), and the DIAMOND table is split accordingly (`split_df`) before parsing.

**Candidate selection (`parseDiamond`).** For each focal protein, subjects are assigned to their species via `strain2species` and the number of distinct species is counted. Two proportions are then computed:

- the proportion of species with an alignment covering the query start (`qstart` < 20 residues);
- the proportion of species with an alignment covering the query end (`qlen − qend` < 20 residues).

A protein is retained as an **N-terminal candidate** when the first proportion is ≤ 5%, and as a **C-terminal candidate** when the second is ≤ 5%. In other words, the terminus under consideration is found in virtually no NR homolog.

> Note: the output column names (`percent_nter_q_15`, `percent_cter_q_15`) suggest a 15-residue threshold. The value actually implemented in `parseDiamond.nf` should be checked.

The module outputs the lists of candidate identifiers (`nter_candidates`, `cter_candidates`) and a per-query statistics table.

### 3.4 Homolog search in the neighboring species (`A3_search_extensions.nf`)

**Protein alignments.** Candidate sequences are extracted from the focal proteome (`seqkit grep`) and aligned against the concatenated local proteome with `ssearch36` (Smith-Waterman), through the `bin/ssearch.sh` wrapper:

```bash
ssearch36 -3 -p -s BL50 -f -11 -g -1 -T${ncpus} -XM${mem}G -m8BCL query subjects
```

i.e. BLOSUM50 matrix, gap open −11 and gap extension −1, BLAST-like tabular output. `blastp` modules (strict and relaxed) are also available as an alternative.

**Homolog filtering (`parseNeighborsBlastp.nf`, Polars).** Alignments with *e*-value ≤ 10⁻⁵, `qcovhsp` ≥ 70% and `ppos` > 50% are retained. Species in which a homolog already has an alignment covering the query terminus ("bad species") are discarded, since no extension is missing there.

**Homolog categorization.** Each query/subject pair is classified according to the position of the alignment on the subject:

| Branch | Variable | "small" category | "big" category |
| --- | --- | --- | --- |
| N-terminal | `sstart` | < 5 | ≥ 5 |
| C-terminal | `slen − send` | < 5 | ≥ 5 |

In the **small** category, the alignment reaches the annotated terminus of the subject: if the missing segment exists, it must lie beyond the annotated CDS, in the *upstream* (N-ter) or *downstream* (C-ter) region of the gene. In the **big** category, the alignment stops before the subject terminus: the segment may be present but divergent within the annotated CDS, or absent.

Results are stored in Parquet format, one table per query, and an `mRNA_species.tsv` table records the transcript → species mapping.

### 3.5 Nucleotide-level analysis of extensions (`modules/nter/`, `modules/cter/`)

#### 3.5.1 "Small" homologs (`small_elongate_and_align`, `bin/{nter,cter}/small_nuc_search.sh`)

**Subject elongation.** For each subject, the elongation length *L* (in nucleotides) is computed from the start (N-ter) or end (C-ter) position of the alignment on the query and from the length of the available intergenic region:

```
L = max(3 × qstart, upstream) × 1.5, rounded down to a multiple of 3
```

The subject CDS coordinates are extracted with `gffread`, then extended by *L* nucleotides with `bedtools slop` on the 5′ side (N-ter) or 3′ side (C-ter), respecting the strand. The extended sequence is retrieved from the indexed concatenated genome.

**Alignments.** Three alignments are performed for each pair:

1. `ssearch36` nucleotide/nucleotide: query CDS against the **elongated** subject CDS;
2. `ssearch36` nucleotide/nucleotide: query CDS against the **non-elongated** subject CDS (control);
3. `ssearch36` protein/protein: query protein against the translation (`faTrans`) of the elongated subject; a `tfasty36` alignment (protein against six-frame translated DNA) is also produced.

The per-subject detection threshold is set to *T* = *L* − 2 nucleotides.

**Parsing (`parse_small_table`).** The elongated alignment is produced in FASTA36 `-m A` format, which gives the position-by-position correspondence between query and subject. From this output the following are computed: the number of gaps on the query and on the subject; the presence of in-frame stop and start codons in the elongated region; the presence of an ATG on the elongated subject facing the query start (`atg_on_elongated_subject_facing_query_start`) and, conversely, of a methionine on the query facing the annotated subject start; the elongation actually recovered, in nucleotides (`raw_recovered_elongation`) and as a fraction of *L* (`recovered_elongation`); and whether the alignment start position on the subject exceeds *L* (`sstart_nuc_gt_elongation`).

#### 3.5.2 "Big" homologs (`ssearch_big_sstart` / `ssearch_big_send`, `bin/{nter,cter}/big_nuc_search.sh`)

Subject CDS are extracted **without elongation** and aligned nucleotide/nucleotide against the query CDS with `ssearch36`. The start (or end) position of the protein alignment is converted to a nucleotide coordinate (×3) and compared with the position obtained from the nucleotide alignment (`parse_big_table`), to determine whether the missing segment is detectable at the DNA level within the annotated CDS.

### 3.6 Summary of parameters and thresholds

| Stage | Parameter | Value |
| --- | --- | --- |
| DIAMOND vs NR | *e*-value; `qcovhsp`; `scovhsp`; `ppos` | ≤ 10⁻⁵; > 60%; > 60%; > 70% |
| DIAMOND vs NR | Taxonomic restriction | Eukaryota (taxid 2759) |
| Candidate selection | Distance to terminus defining a "complete" alignment | < 20 residues (see note, §3.3) |
| Candidate selection | Maximum proportion of species with a complete alignment | ≤ 5% |
| ssearch36 (protein) | Matrix; gap open; gap extension | BLOSUM50; −11; −1 |
| Local homologs | *e*-value; `qcovhsp`; `ppos` | ≤ 10⁻⁵; ≥ 70%; > 50% |
| Categorization | `sstart` or `slen − send` | < 5 (small) / ≥ 5 (big) |
| Elongation | *L* | max(3·`qstart`, upstream) × 1.5, multiple of 3 |
| Elongation | Per-subject threshold *T* | *L* − 2 |
| Parallelization | Number of groups `n` | 30 |

### 3.7 Implementation and reproducibility

The pipeline is written in Nextflow DSL2 and organized into a main workflow (`flows/align.nf`), three subworkflows (`subworkflows/A1_preprocess.nf`, `A2_NR.nf`, `A3_search_extensions.nf`) and atomic modules (`modules/`). Shell scripts in `bin/` encapsulate composite operations (elongation, multiple alignments), and the C++ programs in `src/` (`rmStartStop.cpp`, `parseCoati.cpp`) handle specific sequence processing steps. Intermediate tables are stored in Parquet format and converted to TSV with DuckDB when required. Execution is managed by PBS Pro; AGAT is distributed as a Singularity container (`container/agat.sif`). All executable paths are centralized in `nextflow.config` and exported to the shell scripts as environment variables, with fallback to the `PATH`.

Typical HPC resources:

| Queue | CPUs | RAM | Walltime | Stage |
| --- | --- | --- | --- | --- |
| bim | 70 | 400 GB | 40,000 h | DIAMOND vs NR |
| bim | 16 | 120 GB | 24 h | Alignment parsing |
| bim | 10 | 10 GB | 15 h | ssearch36 |
| common | 4–8 | 8–32 GB | 24 h | Preprocessing |
| lowprio | 4–6 | 4–8 GB | variable | Light parsing |

---

## 4. Outputs

Results are written to `output/nter/` and `output/cter/`:

| File | Content |
| --- | --- |
| `full.blast` | All filtered protein alignments (§3.4) |
| `mRNA_species.tsv` | Transcript → species mapping |
| `bigFinal` | Consolidated results for "big" homologs |
| `smallFinal` | Consolidated results for "small" homologs |

Per-subject intermediate files produced in §3.5.1 are: `*_complete_small_ssearch_elong.tsv` / `.aln` (elongated alignment, tabular and detailed), `*_complete_small_ssearch_short.tsv` (non-elongated alignment), `*_complete_small_tfastx.tsv` (protein alignment) and `*_subjects_thresholds.tsv` (per-subject thresholds *T*).

Columns of `smallFinal` (N-terminal branch):

| Column | Description |
| --- | --- |
| `qseqid` | Query protein identifier (focal species) |
| `sseqid` | Subject homolog identifier |
| `gaps_query` | Number of gaps on the query |
| `gaps_subject` | Number of gaps on the subject |
| `qstart_nuc_elong` | Start position of the elongated nucleotide alignment |
| `sstart_nuc_gt_elongation` | Alignment start on the subject greater than *L* |
| `stop_inframe` | In-frame stop codon in the elongated region |
| `start_inframe` | In-frame start codon in the elongated region |
| `atg_on_elongated_subject_facing_query_start` | ATG on the elongated subject facing the query start |
| `meth_on_query_facing_subject_start` | Methionine on the query facing the annotated subject start |
| `recovered_elongation` | Fraction of *L* recovered by the alignment |
| `raw_recovered_elongation` | Recovered elongation, in nucleotides |
| `category` | "small" or "big" |

---

## 5. Usage

### 5.1 Requirements

```bash
curl -s https://get.nextflow.io | bash
mv nextflow ~/.local/bin/
```

The tools listed in §2.3 must be installed and available in the `PATH`, or their paths set in the configuration.

### 5.2 Configuration (`nextflow.config`)

```groovy
params {
    n             = 30                    // parallelization groups
    output        = "results/"
    focal_dir     = "input/focal/"
    neighbors_dir = "input/neighbors/"
    nr            = "/path/to/nr.dmnd"
    taxdump       = "/path/to/taxdump.tar.gz"
    tmpdir        = "/tmp"

    // Executables (name in PATH or absolute path)
    ssearch = "ssearch36";   tfasty  = "tfasty36"
    seqkit  = "seqkit";      gffread = "gffread"
    faTrans = "faTrans";     faSize  = "faSize"
    diamond = "diamond";     blastp  = "blastp";   makeblastdb = "makeblastdb"
    taxonkit = "taxonkit";   duckdb  = "duckdb"
    samtools = "samtools";   bedtools = "bedtools"
}

process {
    beforeScript = '''
        export SSEARCH="${params.ssearch}"
        export GFFREAD="${params.gffread}"
        # ...
    '''
}
```

Override on the command line or with a configuration file:

```bash
nextflow run main.nf --ssearch /path/to/ssearch36 --diamond /path/to/diamond
nextflow run main.nf -c my_local.config
```

### 5.3 Execution

```bash
# PBS submission
qsub run.sh          # or: qsub nxtflw_run.sh

# Direct execution
nextflow run main.nf -resume
```

### 5.4 Recording software versions

```bash
bash get_versions.sh versions.md
```

---

## 6. Repository layout

```
<pipeline>/
├── main.nf                     # Entry point (align workflow)
├── nextflow.config             # Parameters, executor, tool paths
├── run.sh / nxtflw_run.sh      # PBS submission
├── get_versions.sh             # Software version report
├── flows/align.nf              # Main workflow
├── subworkflows/
│   ├── A1_preprocess.nf        # §3.2
│   ├── A2_NR.nf                # §3.3
│   ├── A3_search_extensions.nf # §3.4
│   └── specifics/              # nter/cter × big/small (§3.5)
├── modules/
│   ├── preprocess/             # extract_protein, makedb, faidx, concat, getGeneCoords, check_annotation, ...
│   ├── NR/                     # diamond, parseDiamond, split_ids, split_df
│   ├── nter/  cter/            # treatBig*, treatSmall*
│   ├── ssearch.nf  blastp.nf  grepSeq.nf
│   ├── parseNeighborsBlastp.nf create_taxon_maps.nf  standardGff.nf
├── bin/
│   ├── ssearch.sh              # ssearch36 wrapper
│   ├── gff.sh  split_dataframe.sh
│   ├── nter/{big,small}_nuc_search.sh
│   └── cter/{big,small}_nuc_search.sh
├── src/                        # rmStartStop.cpp, parseCoati.cpp, json.hpp, compil.sh
├── container/agat.sif
├── input/
│   ├── focal/                  # [Focal].{fna,gff}
│   ├── neighbors/              # [Species].{fna,gff}
│   └── neighbors.genome        # Index of neighbor genomes
├── output/{nter,cter}/
└── test/
```

---

## 7. References

- Buchfink B, Reuter K, Drost HG. Sensitive protein alignments at tree-of-life scale using DIAMOND. *Nat Methods*. 2021;18:366–368.
- Camacho C, Coulouris G, Avagyan V, et al. BLAST+: architecture and applications. *BMC Bioinformatics*. 2009;10:421.
- Cock PJA, Antao T, Chang JT, et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*. 2009;25:1422–1423.
- Dainat J. AGAT: Another Gff Analysis Toolkit to handle annotations in any GTF/GFF format. Zenodo. doi:10.5281/zenodo.3552717.
- Danecek P, Bonfield JK, Liddle J, et al. Twelve years of SAMtools and BCFtools. *GigaScience*. 2021;10:giab008.
- Di Tommaso P, Chatzou M, Floden EW, et al. Nextflow enables reproducible computational workflows. *Nat Biotechnol*. 2017;35:316–319.
- Kent WJ, Sugnet CW, Furey TS, et al. The human genome browser at UCSC. *Genome Res*. 2002;12:996–1006.
- Kurtzer GM, Sochat V, Bauer MW. Singularity: scientific containers for mobility of compute. *PLoS ONE*. 2017;12:e0177459.
- Pearson WR, Lipman DJ. Improved tools for biological sequence comparison. *Proc Natl Acad Sci USA*. 1988;85:2444–2448.
- Pearson WR. Searching protein sequence libraries: comparison of the sensitivity and selectivity of the Smith-Waterman and FASTA algorithms. *Genomics*. 1991;11:635–650.
- Pertea G, Pertea M. GFF Utilities: GffRead and GffCompare. *F1000Research*. 2020;9:304.
- Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing genomic features. *Bioinformatics*. 2010;26:841–842.
- Raasveldt M, Mühleisen H. DuckDB: an embeddable analytical database. *Proc. SIGMOD*. 2019.
- Shen W, Le S, Li Y, Hu F. SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation. *PLoS ONE*. 2016;11:e0163962.
- Shen W, Ren H. TaxonKit: a practical and efficient NCBI taxonomy toolkit. *J Genet Genomics*. 2021;48:844–850.

---

## Changelog

- **v1.0** — Initial release: N- and C-terminal branches, DIAMOND / BLAST+ / FASTA36 integration, PBS Pro parallelization.