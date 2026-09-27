# Runtime Call Graph

Who calls whom, in execution order. This maps the *runtime* path of the data:
`main.nf` → subworkflow → module process → `bin/` script. For the biological
meaning of each step, see `pipeline.md` and `biology.md`; for the shape of the
data passed between steps, see `formats.md`.

Notation: `↳` = "invokes". Files are given as `path :: process` where useful.

## Legend of layers

- **main.nf** — top-level orchestration, calls 4 subworkflows in order.
- **subworkflows/** — group processes into a phase, wire channels together.
- **modules/** — one Nextflow `process` each; a process shells out to a `bin/` script.
- **bin/** — the actual Python/Bash that transforms the data.

---

## 0. Entry point — main.nf

```
main.nf
 ├─ preprocess_input_data()      (A1_preprocess.nf)   → gff, fna, index, proteomes, local_db, geneCoords
 ├─ create_taxon_maps()          (modules/create_taxon_maps.nf) → strain2species, eukaryotes
 ├─ NR(focal_proteome, …)        (A2_NR.nf)           → nter_candidates_IDs, cter_candidates_IDs
 └─ search_homologs(…)           (A3_search_extensions.nf)
```

---

## 1. preprocess_input_data — subworkflows/A1_preprocess.nf

Reads `input/focal/*` and `input/neighbors/*`, produces the shared reference data.

```
A1_preprocess.nf
 ├─ concat_gff / concat_fasta   (preprocess/concat.nf)         ↳ concatenate GFFs / FNAs
 ├─ faidx                       (preprocess/faidx.nf)          ↳ samtools faidx  → index
 ├─ getGeneCoords               (preprocess/getGeneCoords.nf)  ↳ bin/get_gene_coords.py → geneCoords
 ├─ extract_focal_proteome      (preprocess/extract_protein.nf)↳ gffread → focal .faa
 ├─ extract_neighbors_proteome  (preprocess/extract_protein.nf)↳ gffread → neighbor .faa
 ├─ concat_prot                 (preprocess/concat.nf)         ↳ concat neighbor proteins
 └─ blast_makedb                (preprocess/makedb.nf)         ↳ local_db
```

Emits: `full_gff, full_fna, index, focal_proteome, local_db, geneCoords, local_proteome`.

---

## 2. NR — subworkflows/A2_NR.nf

Diamond vs the NR database, then classify each focal protein as N-ter / C-ter candidate.

```
A2_NR.nf
 ├─ diamond (big_diamond)   (NR/diamond.nf)       ↳ diamond blastp focal vs NR → df
 ├─ split_ids              (NR/split_ids.nf)      ↳ chunk ids for parallelism
 ├─ split_df               (NR/split_df.nf)       ↳ split df by chunk
 └─ parseDiamond           (NR/parseDiamond.nf)   ↳ bin/parse_diamond.py
                                                     → per-chunk nter / cter / stats
      collectFile → nter_candidates_IDs, cter_candidates_IDs   (stored in test/)
```

Emits: `nter_candidates_IDs, cter_candidates_IDs`.

---

## 3. search_homologs — subworkflows/A3_search_extensions.nf

Two symmetric branches. **N-ter runs; C-ter is currently commented out**
(A3_search_extensions.nf:61-66) — so nothing under `bin/cter/` executes in the
pipeline as it stands.

### 3a. N-ter branch (active)

```
A3_search_extensions.nf
 ├─ grep_nter          (grepSeq.nf :: seqkitGrep)  ↳ extract candidate seqs from focal proteome
 ├─ ssearch_nter       (blastp.nf :: ssearch)      ↳ bin/ssearch.sh  candidates vs local proteome
 ├─ parseNter          (parseNeighborsBlastp.nf :: nter)  ↳ bin/parse_blastp_nter.py
 │                        → small_sstart (parquet), big_sstart (parquet), bad_species, full_blast
 │
 ├─ nter_small_sstart(small_sstart, …)   (specifics/nter_small_sstart.nf)
 │     ├─ small_elongate_and_align   (nter/treatSmallSstart.nf :: small_elongate_and_align)
 │     │                               ↳ bin/nter/small_nuc_search.sh
 │     │                                 → *_subjects_thresholds.tsv, ssearch tables, sequences
 │     └─ parse_small_table          (nter/treatSmallSstart.nf :: parse_small_table)
 │                                     ↳ bin/nter/parse_small_table.py
 │                                       → {qseqid}_final_res.tsv
 │           collectFile → output/nter/smallFinal
 │
 └─ nter_big_sstart(big_sstart, …)       (specifics/nter_big_sstart.nf)
       ├─ ssearch_big_sstart         (nter/treatBigSstart.nf)  ↳ bin/nter/big_nuc_search.sh
       │                                                          → ssearch tsv + qstarts.tsv
       └─ parse_big_table            (nter/treatBigSstart.nf)  ↳ bin/nter/parse_big_table.py
                                                                 → {qseqid}_final_res.tsv
             collectFile → output/nter/bigFinal
```

### 3b. C-ter branch (commented out — mirror of N-ter)

Same shape, `send/slend` instead of `sstart`, downstream instead of upstream.
To enable, uncomment A3_search_extensions.nf:61-66.

```
 ├─ grep_cter          (grepSeq.nf :: seqkitGrep)
 ├─ ssearch_cter       (blastp.nf :: ssearch)
 ├─ parseCter          (parseNeighborsBlastp.nf :: cter)  ↳ bin/parse_blastp_cter.py
 │                        → small_send, big_send, bad_species, full_blast
 │
 ├─ cter_small_send(small_send, …)   (specifics/cter_small_send.nf)
 │     ├─ small_elongate_and_align   (cter/treatSmallSend.nf)  ↳ bin/cter/small_nuc_search.sh
 │     │                                                          → *_subjects_thresholds.tsv (now
 │     │                                                            carries elongated_length)
 │     └─ parse_small_table          (cter/treatSmallSend.nf)  ↳ bin/cter/parse_small_table.py
 │                                     [STATUS: incomplete — no per-subject loop yet]
 │           collectFile → output/cter/smallFinal
 │
 └─ cter_big_send(big_send, …)       (specifics/cter_big_send.nf)
       ├─ ssearch_big_send           (cter/treatBigSend.nf)   ↳ bin/cter/big_nuc_search.sh
       │                                                          → ssearch tsv + qends.tsv
       └─ parse_big_table            (cter/treatBigSend.nf)   ↳ bin/cter/parse_big_table.py
             collectFile → output/cter/bigFinal
```

---

## Quick lookup — process → script

| Process (module)                      | bin/ script                          |
| ------------------------------------- | ------------------------------------ |
| getGeneCoords                         | bin/get_gene_coords.py               |
| parseDiamond                          | bin/parse_diamond.py                 |
| ssearch (blastp.nf)                   | bin/ssearch.sh                       |
| nter / cter (parseNeighborsBlastp.nf) | bin/parse_blastp_{nter,cter}.py      |
| small_elongate_and_align (nter/cter)  | bin/{nter,cter}/small_nuc_search.sh  |
| parse_small_table (nter/cter)         | bin/{nter,cter}/parse_small_table.py |
| ssearch_big_* (nter/cter)             | bin/{nter,cter}/big_nuc_search.sh    |
| parse_big_table (nter/cter)           | bin/{nter,cter}/parse_big_table.py   |

## Final outputs

- output/nter/{smallFinal,bigFinal}  ← active
- output/cter/{smallFinal,bigFinal}  ← produced only once branch 3b is uncommented
- output/{nter,cter}/mRNA_species.tsv ← sseqid → species map (needed to attach
  neighbor results to tree leaves; see the tree/triage discussion)
