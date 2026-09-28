# DMEL — Détection et validation d'extensions terminales protéiques chez *Drosophila melanogaster*

Pipeline Nextflow pour l'identification, chez *D. melanogaster* (Dmel), de segments N-terminaux et C-terminaux absents des homologues annotés des autres espèces, et pour la recherche de ces segments dans l'ADN génomique non annoté des drosophiles voisines.

*Pipeline développé par Simon Herman — équipe BIM.*

## Sommaire

1. [Contexte et objectif](#1-contexte-et-objectif)
2. [Matériel](#2-matériel)
3. [Méthodes](#3-méthodes)
4. [Sorties](#4-sorties)
5. [Utilisation](#5-utilisation)
6. [Organisation du dépôt](#6-organisation-du-dépôt)
7. [Références](#7-références)

---

## 1. Contexte et objectif

L'annotation automatique des génomes produit fréquemment des codons d'initiation ou de terminaison erronés, ce qui conduit à ignorer des segments codants pourtant conservés. *D. melanogaster* étant l'espèce du genre la mieux annotée, un certain nombre de ses protéines apparaissent plus longues que leurs homologues annotés chez les autres espèces. Pour chacune de ces différences de longueur, trois hypothèses sont envisageables :

1. une **innovation réelle** de la lignée Dmel ;
2. une **erreur d'annotation de Dmel** (codon START placé trop en amont, ou STOP trop en aval) ;
3. une **sous-annotation de l'espèce voisine**, le segment étant présent dans son génome mais non prédit comme codant.

Le pipeline vise à discriminer ces hypothèses en trois temps : (i) identification, par recherche contre la base NR, des protéines de Dmel dont une extrémité n'est retrouvée chez pratiquement aucun homologue ; (ii) alignement de ces candidats contre les protéomes des drosophiles voisines ; (iii) pour chaque paire candidat/homologue, extension virtuelle du gène voisin dans l'ADN génomique et alignement nucléotidique afin de tester la présence du segment manquant.

---

## 2. Matériel

### 2.1 Espèce focale

| Fichier | Contenu | Taille |
| --- | --- | --- |
| `input/focal/Dmel.fna` | Séquences génomiques de *D. melanogaster* (assemblage : *à préciser*) | ~145 Mo |
| `input/focal/Dmel.gff` | Annotation GFF3 associée (version : *à préciser*) | ~164 Mo |
| `input/focal/Dmel.proteome` | Protéome dérivé de l'annotation | ~20 Mo |

### 2.2 Espèces voisines

Soixante-trois espèces du genre *Drosophila* (`input/neighbors/`), chacune représentée par un génome (`[Species].fna`) et son annotation GFF3 (`[Species].gff`). L'index `input/neighbors.genome` recense les génomes disponibles. Les relations phylogénétiques entre espèces sont décrites dans `Drosophila_tree.nwk` (format Newick) et incluent notamment les espèces du groupe *melanogaster* (Dsim, Dsec, Dmau, Dere, Dsan, Dyak, Dtei, Deug).

### 2.3 Bases de données externes

- **NCBI NR** au format DIAMOND (`nr_2.0.13.dmnd`), interrogée avec une restriction aux eucaryotes (taxid 2759).
- **NCBI Taxonomy** (`taxdump`), utilisée pour la conversion souche → espèce et pour la définition de l'ensemble des taxids eucaryotes.

### 2.4 Outils logiciels

Tous les chemins d'exécutables sont paramétrables dans `nextflow.config` (voir §5.2). Les versions exactes utilisées sont à renseigner dans la colonne correspondante.

| Outil | Version | Usage dans le pipeline | Référence |
| --- | --- | --- | --- |
| Nextflow (DSL2) | — | Orchestration du workflow | Di Tommaso et al., 2017 |
| Singularity | — | Conteneurisation (AGAT) | Kurtzer et al., 2017 |
| PBS Pro | — | Ordonnancement HPC | — |
| DIAMOND | — | Recherche de similarité contre NR | Buchfink et al., 2021 |
| BLAST+ (`makeblastdb`, `blastp`) | — | Base de données protéique locale, alignements protéiques | Camacho et al., 2009 |
| FASTA36 (`ssearch36`, `tfasty36`) | — | Alignements Smith-Waterman prot/prot, nuc/nuc et prot/ADN traduit | Pearson & Lipman, 1988 ; Pearson, 1991 |
| gffread | — | Extraction des CDS et des protéines depuis GFF + génome | Pertea & Pertea, 2020 |
| BEDTools (`slop`) | — | Extension de coordonnées génomiques | Quinlan & Hall, 2010 |
| SAMtools (`faidx`) | — | Indexation des génomes | Danecek et al., 2021 |
| SeqKit | — | Manipulation de fichiers FASTA | Shen et al., 2016 |
| TaxonKit | — | Manipulation de la taxonomie NCBI | Shen & Ren, 2021 |
| AGAT | — | Standardisation des GFF | Dainat |
| UCSC utilities (`faTrans`, `faSize`) | — | Traduction ADN → protéine, tailles de séquences | Kent et al., 2002 |
| DuckDB | — | Conversion Parquet → TSV | Raasveldt & Mühleisen, 2019 |
| Python (Polars, Biopython, gff3_parser) | — | Parsing et traitement des tables | Cock et al., 2009 |

---

## 3. Méthodes

### 3.1 Vue d'ensemble

Le pipeline s'articule en quatre étapes (Figure 1) : prétraitement des données (§3.2), identification des candidats par recherche contre NR (§3.3), recherche d'homologues chez les espèces voisines (§3.4) et analyse nucléotidique des extensions (§3.5). Les extrémités N-terminale (« nter ») et C-terminale (« cter ») sont traitées par deux branches symétriques.

```
                              main.nf → align
                                    │
             ┌──────────────────────┼──────────────────────┐
             ▼                      ▼                      ▼
     A1_preprocess          create_taxon_maps            A2_NR
  extraction protéomes      strain2species        DIAMOND vs NR (euk.)
  concat / faidx / makedb   eukaryotes taxids     parseDiamond → candidats
  getGeneCoords                                   nter / cter
             │                                              │
             └──────────────────────┬───────────────────────┘
                                    ▼
                         A3_search_extensions
                  ssearch36 candidats vs protéome local
                  parseNter / parseCter → big | small
                                    │
                     ┌──────────────┴──────────────┐
                     ▼                             ▼
             big_sstart / big_send       small_sstart / small_send
             ssearch nuc (CDS seul)      élongation (bedtools slop)
             + parse                     + alignements multiples
                                         + parse (gaps, codons, élongation récupérée)
                                    │
                                    ▼
                          output/nter/  output/cter/
```

*Figure 1. Architecture du pipeline.*

### 3.2 Prétraitement des données (`A1_preprocess.nf`, `create_taxon_maps.nf`)

**Extraction des protéomes.** Pour chaque espèce, le protéome est dérivé de l'annotation et du génome avec `gffread` (option `-J`, qui écarte les transcrits sans codon START ou avec un STOP en phase), puis les gaps sont retirés avec `seqkit` :

```bash
gffread -J -y - -g ${fna} ${gff} | seqkit seq --remove-gaps - > ${species}.faa
```

**Constitution des jeux de données concaténés.** Les GFF, génomes et protéomes de l'ensemble des espèces voisines sont concaténés (`full_gff`, `full_fna`, `local_proteome`). Le génome concaténé est indexé avec `samtools faidx`, et une base de données protéique locale est construite avec `makeblastdb` (`-dbtype prot -parse_seqids -hash_index`).

**Coordonnées géniques.** Pour chaque gène, le module `getGeneCoords` calcule, en tenant compte du brin, la distance entre le début du CDS et le début du gène (région *upstream*) et entre la fin du CDS et la fin du gène (région *downstream*). Ces distances bornent l'élongation appliquée en §3.5.1.

**Tables taxonomiques.** À partir du `taxdump` NCBI, `taxonkit` produit une table de correspondance souche → espèce (`strain2species.csv`) et la liste des taxids eucaryotes descendant du taxid 2759 (`eukaryotes.csv`).

### 3.3 Identification des candidats par recherche contre NR (`A2_NR.nf`)

**Recherche de similarité.** Le protéome de Dmel est aligné contre NR avec `diamond blastp` en mode `--fast`, restreint aux eucaryotes, sans limite sur le nombre de cibles :

```bash
diamond blastp --query ${query} --db ${db} --taxonlist 2759 \
  --outfmt 6 qseqid sseqid qlen qstart qend qcovhsp scovhsp ppos staxids \
  --max-target-seqs 0 --fast -e 0.00001
```

Les alignements sont ensuite filtrés sur la couverture de la query (`qcovhsp` > 60 %), la couverture du sujet (`scovhsp` > 60 %) et le pourcentage de positions positives (`ppos` > 70 %).

**Parallélisation.** Les identifiants de queries sont répartis en `n` groupes (30 par défaut ; `split_ids`), et la table DIAMOND est découpée en conséquence (`split_df`) avant parsing.

**Sélection des candidats (`parseDiamond`).** Pour chaque protéine de Dmel, les sujets sont rattachés à leur espèce via `strain2species` et le nombre d'espèces distinctes est comptabilisé. Deux proportions sont alors calculées :

- la proportion d'espèces possédant un alignement couvrant le début de la query (`qstart` < 20 résidus) ;
- la proportion d'espèces possédant un alignement couvrant la fin de la query (`qlen − qend` < 20 résidus).

Une protéine est retenue comme **candidat N-terminal** lorsque la première proportion est ≤ 5 %, et comme **candidat C-terminal** lorsque la seconde est ≤ 5 %. Autrement dit, l'extrémité considérée de la protéine de Dmel n'est retrouvée chez pratiquement aucun homologue de NR.

> Note : le nom des colonnes de sortie (`percent_nter_q_15`, `percent_cter_q_15`) suggère un seuil de 15 résidus. La valeur effectivement codée dans `parseDiamond.nf` est à vérifier.

Le module produit la liste des identifiants candidats (`nter_candidates`, `cter_candidates`) ainsi qu'une table de statistiques par query.

### 3.4 Recherche d'homologues chez les espèces voisines (`A3_search_extensions.nf`)

**Alignements protéiques.** Les séquences des candidats sont extraites du protéome focal (`seqkit grep`) puis alignées contre le protéome local concaténé avec `ssearch36` (Smith-Waterman), via le wrapper `bin/ssearch.sh` :

```bash
ssearch36 -3 -p -s BL50 -f -11 -g -1 -T${ncpus} -XM${mem}G -m8BCL query subjects
```

soit une matrice BLOSUM50, des pénalités de gap d'ouverture −11 et d'extension −1, et une sortie tabulaire de type BLAST. Des modules `blastp` (strict et « souple ») sont également disponibles en alternative.

**Filtrage des homologues (`parseNeighborsBlastp.nf`, Polars).** Sont conservés les alignements satisfaisant *e*-value ≤ 10⁻⁵, `qcovhsp` ≥ 70 % et `ppos` > 50 %. Les espèces dont un homologue présente déjà un alignement couvrant l'extrémité de la query (« bad species ») sont écartées, l'extension n'y étant pas manquante.

**Catégorisation des homologues.** Chaque paire query/sujet est classée selon la position de l'alignement sur le sujet :

| Branche | Variable | Catégorie « small » | Catégorie « big » |
| --- | --- | --- | --- |
| N-terminale | `sstart` | < 5 | ≥ 5 |
| C-terminale | `slen − send` | < 5 | ≥ 5 |

Dans la catégorie **small**, l'alignement atteint l'extrémité annotée du sujet : si le segment manquant existe, il doit se situer au-delà du CDS annoté, dans la région *upstream* (N-ter) ou *downstream* (C-ter) du gène. Dans la catégorie **big**, l'alignement s'interrompt avant l'extrémité du sujet : le segment peut être présent mais divergent au sein du CDS annoté, ou absent.

Les résultats sont stockés au format Parquet, une table par query, et une table `mRNA_species.tsv` conserve la correspondance transcrit → espèce.

### 3.5 Analyse nucléotidique des extensions (`modules/nter/`, `modules/cter/`)

#### 3.5.1 Homologues de catégorie « small » (`small_elongate_and_align`, `bin/{nter,cter}/small_nuc_search.sh`)

**Élongation des sujets.** Pour chaque sujet, la longueur d'élongation *L* (en nucléotides) est calculée à partir de la position de début (N-ter) ou de fin (C-ter) de l'alignement sur la query et de la longueur de la région intergénique disponible :

```
L = max(3 × qstart, upstream) × 1.5, arrondi au multiple de 3 inférieur
```

Les coordonnées du CDS du sujet sont extraites avec `gffread`, puis étendues de *L* nucléotides avec `bedtools slop` du côté 5′ (N-ter) ou 3′ (C-ter), en respectant le brin. La séquence étendue est extraite du génome concaténé indexé.

**Alignements.** Trois alignements sont réalisés pour chaque paire :

1. `ssearch36` nucléotide/nucléotide : CDS de la query contre le CDS du sujet **élongé** ;
2. `ssearch36` nucléotide/nucléotide : CDS de la query contre le CDS du sujet **non élongé** (contrôle) ;
3. `ssearch36` protéine/protéine : protéine de la query contre la traduction (`faTrans`) du sujet élongé ; un alignement `tfasty36` (protéine contre ADN traduit dans les six cadres) est également produit.

Le seuil de détection propre à chaque sujet est fixé à *T* = *L* − 2 nucléotides.

**Parsing (`parse_small_table`).** L'alignement élongé est produit au format `-m A` de FASTA36, qui fournit la correspondance position par position entre query et sujet. À partir de cette sortie sont calculés : le nombre de gaps sur la query et sur le sujet ; la présence de codons STOP et START en phase dans la région élongée ; la présence d'un ATG sur le sujet élongé en regard du début de la query (`atg_on_elongated_subject_facing_query_start`) et, réciproquement, d'une méthionine sur la query en regard du début du sujet annoté ; l'élongation effectivement récupérée, en nucléotides (`raw_recovered_elongation`) et rapportée à *L* (`recovered_elongation`) ; et le dépassement éventuel de *L* par la position de début d'alignement sur le sujet (`sstart_nuc_gt_elongation`).

#### 3.5.2 Homologues de catégorie « big » (`ssearch_big_sstart` / `ssearch_big_send`, `bin/{nter,cter}/big_nuc_search.sh`)

Les CDS des sujets sont extraits **sans élongation**, puis alignés en nucléotide/nucléotide contre le CDS de la query avec `ssearch36`. La position de début (ou de fin) de l'alignement protéique est convertie en coordonnée nucléotidique (×3) et comparée à la position obtenue par l'alignement nucléotidique (`parse_big_table`), afin de déterminer si le segment manquant est détectable au niveau ADN au sein du CDS annoté.

### 3.6 Récapitulatif des paramètres et seuils

| Étape | Paramètre | Valeur |
| --- | --- | --- |
| DIAMOND vs NR | *e*-value ; `qcovhsp` ; `scovhsp` ; `ppos` | ≤ 10⁻⁵ ; > 60 % ; > 60 % ; > 70 % |
| DIAMOND vs NR | Restriction taxonomique | Eucaryotes (taxid 2759) |
| Sélection des candidats | Distance à l'extrémité définissant un alignement « complet » | < 20 résidus (voir note §3.3) |
| Sélection des candidats | Proportion maximale d'espèces avec alignement complet | ≤ 5 % |
| ssearch36 (prot/prot) | Matrice ; gap ouverture ; gap extension | BLOSUM50 ; −11 ; −1 |
| Homologues locaux | *e*-value ; `qcovhsp` ; `ppos` | ≤ 10⁻⁵ ; ≥ 70 % ; > 50 % |
| Catégorisation | `sstart` ou `slen − send` | < 5 (small) / ≥ 5 (big) |
| Élongation | *L* | max(3·`qstart`, upstream) × 1,5, multiple de 3 |
| Élongation | Seuil par sujet *T* | *L* − 2 |
| Parallélisation | Nombre de groupes `n` | 30 |

### 3.7 Implémentation et reproductibilité

Le pipeline est écrit en Nextflow DSL2 et organisé en un workflow principal (`flows/align.nf`), trois sous-workflows (`subworkflows/A1_preprocess.nf`, `A2_NR.nf`, `A3_search_extensions.nf`) et des modules atomiques (`modules/`). Les scripts shell de `bin/` encapsulent les opérations composites (élongation, alignements multiples) et les programmes C++ de `src/` (`rmStartStop.cpp`, `parseCoati.cpp`) assurent des traitements ponctuels sur les séquences. Les tables intermédiaires sont stockées au format Parquet et converties en TSV avec DuckDB lorsque nécessaire. L'exécution est prise en charge par PBS Pro ; AGAT est distribué sous forme de conteneur Singularity (`container/agat.sif`). Les chemins de tous les exécutables sont centralisés dans `nextflow.config` et exportés vers les scripts shell sous forme de variables d'environnement, avec repli sur le `PATH`.

Ressources HPC typiques :

| Queue | CPUs | RAM | Walltime | Étape |
| --- | --- | --- | --- | --- |
| bim | 70 | 400 Go | 40 000 h | DIAMOND vs NR |
| bim | 16 | 120 Go | 24 h | Parsing des alignements |
| bim | 10 | 10 Go | 15 h | ssearch36 |
| common | 4–8 | 8–32 Go | 24 h | Prétraitement |
| lowprio | 4–6 | 4–8 Go | variable | Parsing léger |

---

## 4. Sorties

Les résultats sont écrits dans `output/nter/` et `output/cter/` :

| Fichier | Contenu |
| --- | --- |
| `full.blast` | Ensemble des alignements protéiques filtrés (§3.4) |
| `mRNA_species.tsv` | Correspondance transcrit → espèce |
| `bigFinal` | Résultats consolidés des homologues de catégorie « big » |
| `smallFinal` | Résultats consolidés des homologues de catégorie « small » |

Les fichiers intermédiaires produits par sujet en §3.5.1 sont : `*_complete_small_ssearch_elong.tsv` / `.aln` (alignement élongé, tabulaire et détaillé), `*_complete_small_ssearch_short.tsv` (alignement non élongé), `*_complete_small_tfastx.tsv` (alignement protéique) et `*_subjects_thresholds.tsv` (seuils *T* par sujet).

Colonnes de `smallFinal` (branche N-terminale) :

| Colonne | Description |
| --- | --- |
| `qseqid` | Identifiant de la protéine query (Dmel) |
| `sseqid` | Identifiant de l'homologue sujet |
| `gaps_query` | Nombre de gaps sur la query |
| `gaps_subject` | Nombre de gaps sur le sujet |
| `qstart_nuc_elong` | Position de début de l'alignement nucléotidique élongé |
| `sstart_nuc_gt_elongation` | Début d'alignement sur le sujet supérieur à *L* |
| `stop_inframe` | Codon STOP en phase dans la région élongée |
| `start_inframe` | Codon START en phase dans la région élongée |
| `atg_on_elongated_subject_facing_query_start` | ATG sur le sujet élongé en regard du début de la query |
| `meth_on_query_facing_subject_start` | Méthionine sur la query en regard du début du sujet annoté |
| `recovered_elongation` | Fraction de *L* récupérée par l'alignement |
| `raw_recovered_elongation` | Élongation récupérée, en nucléotides |
| `category` | « small » ou « big » |

---

## 5. Utilisation

### 5.1 Prérequis

```bash
curl -s https://get.nextflow.io | bash
mv nextflow ~/.local/bin/
```

Les outils listés en §2.4 doivent être installés et accessibles dans le `PATH`, ou leurs chemins renseignés dans la configuration.

### 5.2 Configuration (`nextflow.config`)

```groovy
params {
    n             = 30                    // groupes de parallélisation
    output        = "results/"
    focal_dir     = "input/focal/"
    neighbors_dir = "input/neighbors/"
    nr            = "/path/to/nr.dmnd"
    taxdump       = "/path/to/taxdump.tar.gz"
    tmpdir        = "/tmp"

    // Exécutables (nom dans le PATH ou chemin absolu)
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

Surcharge en ligne de commande ou par fichier de configuration :

```bash
nextflow run main.nf --ssearch /path/to/ssearch36 --diamond /path/to/diamond
nextflow run main.nf -c my_local.config
```

### 5.3 Exécution

```bash
# Soumission PBS
qsub run.sh          # ou : qsub nxtflw_run.sh

# Exécution directe
nextflow run main.nf -resume
```

---

## 6. Organisation du dépôt

```
DMEL/
├── main.nf                     # Point d'entrée (workflow align)
├── nextflow.config             # Paramètres, exécuteur, chemins d'outils
├── run.sh / nxtflw_run.sh      # Soumission PBS
├── flows/align.nf              # Workflow principal
├── subworkflows/
│   ├── A1_preprocess.nf        # §3.2
│   ├── A2_NR.nf                # §3.3
│   ├── A3_search_extensions.nf # §3.4
│   └── specifics/              # nter/cter × big/small (§3.5)
├── modules/
│   ├── preprocess/             # extract_protein, makedb, faidx, concat, getGeneCoords, ...
│   ├── NR/                     # diamond, parseDiamond, split_ids, split_df
│   ├── nter/  cter/            # treatBig*, treatSmall*
│   ├── ssearch.nf  blastp.nf  grepSeq.nf
│   ├── parseNeighborsBlastp.nf create_taxon_maps.nf  standardGff.nf
├── bin/
│   ├── ssearch.sh              # Wrapper ssearch36
│   ├── gff.sh  split_dataframe.sh
│   ├── nter/{big,small}_nuc_search.sh
│   └── cter/{big,small}_nuc_search.sh
├── src/                        # rmStartStop.cpp, parseCoati.cpp, json.hpp, compil.sh
├── container/agat.sif
├── input/
│   ├── focal/                  # Dmel.{fna,gff,proteome}
│   ├── neighbors/              # [Species].{fna,gff}
│   └── neighbors.genome
├── output/{nter,cter}/
├── test/
├── Drosophila_tree.nwk
└── Dmel.tsv                    # Table DIAMOND vs NR (~1,1 Go)
```

---

## 7. Références

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

## Historique

- **v1.0** — Pipeline initial pour *D. melanogaster* : branches N- et C-terminales, intégration DIAMOND / BLAST+ / FASTA36, parallélisation PBS Pro.