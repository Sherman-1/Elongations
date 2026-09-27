# DMEL - Pipeline de Détection d'Extensions de Protéines chez *Drosophila melanogaster*

## 📋 Table des Matières

1. [Vue d'ensemble](#vue-densemble)
2. [Objectif Scientifique](#objectif-scientifique)
3. [Architecture du Pipeline](#architecture-du-pipeline)
4. [Structure du Projet](#structure-du-projet)
5. [Données d'Entrée](#données-dentrée)
6. [Workflow Détaillé](#workflow-détaillé)
7. [Modules Nextflow](#modules-nextflow)
8. [Scripts Shell (bin/)](#scripts-shell-bin)
9. [Outils Utilisés](#outils-utilisés)
10. [Configuration et Exécution](#configuration-et-exécution)
11. [Outputs](#outputs)
12. [Reconstruction du Pipeline](#reconstruction-du-pipeline)

---

## Vue d'ensemble

Ce pipeline Nextflow est conçu pour **valider et caractériser les extensions** N-terminales (nter) et C-terminales (cter) des protéines de *Drosophila melanogaster* (Dmel). Pour cela, il utilise Dmel comme référence de haute qualité et recherche si ses extensions "uniques" sont conservées dans l'ADN génomique (non annoté) des espèces voisines.

**Technologies principales :**
- **Nextflow** : Orchestration du workflow 
- **Singularity** : Conteneurisation
- **PBS Pro** : Gestionnaire de ressources HPC
- **Python/Polars** : Traitement de données haute performance
- **FASTA36 (ssearch36/tfastx)** : Alignements Smith-Waterman précis

---

## Objectif Scientifique

### Problématique
L'annotation automatique des génomes conduit souvent à des codons START ou STOP erronés, ignorant parfois des segments codants conservés. *D. melanogaster*, étant l'espèce la mieux annotée, possède des protéines qui semblent plus longues que leurs homologues chez d'autres espèces. Ce pipeline cherche à déterminer si cette "longueur supplémentaire" est :
- Une **innovation réelle** de Dmel.
- Une **erreur d'annotation de Dmel** (START trop en amont).
- Ou, plus fréquemment, une **sous-annotation des autres espèces** (le segment existe dans leur ADN mais n'a pas été prédit comme codant).

### Approche
1. **Recherche contre NR (Non-Redundant database)** : Identifier les protéines de Dmel qui possèdent un segment (N-ter ou C-ter) non retrouvé chez 95% des autres espèces.
2. **Recherche locale ciblée** : Comparer ces segments de Dmel contre les génomes des drosophiles voisines.
3. **Analyse nucléotidique (Slop)** : Étendre virtuellement les gènes des voisins (sujets) et aligner la séquence de Dmel (query) pour voir si le segment y est présent.

### Critères de sélection des candidats (via Diamond vs NR)
- **N-ter** : Protéines de Dmel où **moins de 5%** des homologues dans NR possèdent un alignement couvrant le début de la protéine (qstart < 20). En d'autres termes, Dmel a un début "unique".
- **C-ter** : Protéines de Dmel où **moins de 5%** des homologues dans NR possèdent un alignement couvrant la fin de la protéine (qlend < 20).

---

## Architecture du Pipeline

```
┌─────────────────────────────────────────────────────────────────┐
│                          main.nf                                 │
│                     └── align workflow                          │
└─────────────────────────────────────────────────────────────────┘
                                │
        ┌───────────────────────┼───────────────────────┐
        ▼                       ▼                       ▼
┌───────────────┐    ┌──────────────────┐    ┌──────────────────┐
│ A1_preprocess │    │  create_taxon_   │    │     A2_NR        │
│               │    │     maps         │    │                  │
│ - extract     │    │ - strain2species │    │ - diamond vs NR  │
│   proteomes   │    │ - eukaryotes     │    │ - parse results  │
│ - concat      │    │   taxids         │    │ - identify       │
│ - makedb      │    │                  │    │   candidates     │
│ - faidx       │    └──────────────────┘    └────────┬─────────┘
│ - geneCoords  │                                      │
└───────┬───────┘                                      │
        │                                              │
        │          ┌───────────────────────────────────┘
        │          │
        ▼          ▼
┌─────────────────────────────────────────────────────────────────┐
│                    A3_search_extensions                          │
│                                                                  │
│  ┌─────────────────────────┐    ┌─────────────────────────────┐  │
│  │     N-terminal          │    │       C-terminal            │  │
│  │                         │    │                             │  │
│  │  grep_nter              │    │    grep_cter                │  │
│  │  ssearch_nter           │    │    ssearch_cter             │  │
│  │  parseNter              │    │    parseCter                │  │
│  │         │               │    │           │                 │  │
│  │    ┌────┴────┐          │    │      ┌────┴────┐            │  │
│  │    ▼         ▼          │    │      ▼         ▼            │  │
│  │ big_sstart small_sstart │    │  big_send   small_send      │  │
│  │    │         │          │    │      │         │            │  │
│  │    ▼         ▼          │    │      ▼         ▼            │  │
│  │ ssearch   elongate_     │    │  ssearch   elongate_        │  │
│  │ + parse   + align       │    │  + parse   + align          │  │
│  │           + parse       │    │            + parse          │  │
│  └─────────────────────────┘    └─────────────────────────────┘  │
└─────────────────────────────────────────────────────────────────┘
                                │
                                ▼
                         ┌──────────────┐
                         │   OUTPUTS    │
                         │              │
                         │ output/nter/ │
                         │ output/cter/ │
                         └──────────────┘
```

---

## Structure du Projet

```
DMEL/
├── main.nf                    # Point d'entrée du workflow
├── nextflow.config            # Configuration Nextflow (params, executor, singularity)
├── run.sh / nxtflw_run.sh     # Scripts de soumission PBS
│
├── flows/                     # Workflows principaux
│   ├── align.nf               # Workflow principal d'alignement
│   └── test.nf                # Workflow de test
│
├── subworkflows/              # Sous-workflows
│   ├── A1_preprocess.nf       # Prétraitement des données
│   ├── A2_NR.nf               # Recherche contre NR
│   ├── A3_search_extensions.nf # Recherche d'extensions
│   └── specifics/             # Traitement spécifique nter/cter
│       ├── nter_big_sstart.nf
│       ├── nter_small_sstart.nf
│       ├── cter_big_send.nf
│       └── cter_small_send.nf
│
├── modules/                   # Modules Nextflow (processus)
│   ├── preprocess/            # Modules de prétraitement
│   │   ├── extract_protein.nf
│   │   ├── makedb.nf
│   │   ├── faidx.nf
│   │   ├── concat.nf
│   │   ├── sample_proteins.nf
│   │   ├── getGeneCoords.nf
│   │   └── check_annotation.nf
│   ├── NR/                    # Modules pour recherche NR
│   │   ├── diamond.nf
│   │   ├── parseDiamond.nf
│   │   ├── split_ids.nf
│   │   └── split_df.nf
│   ├── nter/                  # Modules N-terminaux
│   │   ├── treatBigSstart.nf
│   │   ├── treatSmallSstart.nf
│   │   └── final.nf
│   ├── cter/                  # Modules C-terminaux
│   │   ├── treatBigSend.nf
│   │   └── treatSmallSend.nf
│   ├── blastp.nf              # Alignements protéiques
│   ├── grepSeq.nf             # Extraction de séquences
│   ├── create_taxon_maps.nf   # Création maps taxonomiques
│   ├── parseNeighborsBlastp.nf # Parsing des résultats blastp
│   ├── standardGff.nf         # Standardisation GFF
│   └── ssearch.nf             # Alignements ssearch
│
├── bin/                       # Scripts exécutables
│   ├── ssearch.sh             # Wrapper ssearch36
│   ├── gff.sh                 # Utilitaires GFF
│   ├── split_dataframe.sh     # Split dataframes
│   ├── test.sh                # Tests
│   ├── nter/                  # Scripts N-terminaux
│   │   ├── big_nuc_search.sh
│   │   └── small_nuc_search.sh
│   └── cter/                  # Scripts C-terminaux
│       ├── big_nuc_search.sh
│       └── small_nuc_search.sh
│
├── src/                       # Code source C++
│   ├── rmStartStop.cpp        # Suppression start/stop codons
│   ├── parseCoati.cpp         # Parsing alignements Coati
│   ├── json.hpp               # Bibliothèque JSON
│   └── compil.sh              # Script de compilation
│
├── container/                 # Conteneurs Singularity
│   └── agat.sif               # Container AGAT
│
├── input/                     # Données d'entrée
│   ├── focal/                 # Espèce focale (D. melanogaster)
│   │   ├── Dmel.fna           # Séquences génomiques
│   │   ├── Dmel.gff           # Annotations
│   │   └── Dmel.proteome      # Protéome
│   ├── neighbors/             # Espèces voisines (~63 espèces)
│   │   └── [Species].{fna,gff}
│   └── neighbors.genome       # Index des génomes voisins
│
├── output/                    # Résultats
│   ├── nter/                  # Candidats N-terminaux
│   └── cter/                  # Candidats C-terminaux
│
├── test/                      # Fichiers de test et résultats intermédiaires
├── Drosophila_tree.nwk        # Arbre phylogénétique Newick
└── Dmel.tsv                   # Données tabulaires (~1.1GB)
```

---

## Données d'Entrée

### Espèce Focale (input/focal/)
- **Dmel.fna** : Séquences génomiques de *D. melanogaster* (~145 MB)
- **Dmel.gff** : Annotations GFF3 (~164 MB)
- **Dmel.proteome** : Protéome extrait (~20 MB)

### Espèces Voisines (input/neighbors/)
~63 espèces du genre *Drosophila* avec pour chacune :
- `[Species].fna` : Séquences génomiques
- `[Species].gff` : Annotations GFF3

### Arbre Phylogénétique
```newick
(((((Drosophila_rhopaloa:16.32455000,Drosophila_elegans:16.32455000)...
```
Contient les relations phylogénétiques entre les espèces utilisées, incluant :
- Dsim, Dsec, Dmau, Dmel, Dere, Dsan, Dyak, Dtei, Deug...

### Base de Données Externe
- **NR (Non-Redundant)** : `/datas/NR/nr_2.0.13.dmnd` (format Diamond)
- **Taxdump** : Données taxonomiques NCBI

---

## Workflow Détaillé

### Étape 1 : Prétraitement (`A1_preprocess.nf`)

**But** : Préparer les données pour l'analyse

**Processus** :

1. **extract_protein** : Extraction des protéomes depuis GFF+FNA
   ```bash
   gffread -J -y - -g ${fna} ${gff} | seqkit seq --remove-gaps - > ${species}.faa
   ```

2. **concat** : Concaténation de tous les fichiers GFF, FNA, et protéomes

3. **faidx** : Indexation des génomes avec samtools
   ```bash
   samtools faidx ${genome}
   ```

4. **blast_makedb** : Création de la base BLAST locale
   ```bash
   makeblastdb -in $fasta -dbtype 'prot' -out blast -parse_seqids -hash_index
   ```

5. **getGeneCoords** : Calcul des coordonnées upstream/downstream pour chaque gène
   - Calcule la distance entre le début du CDS et le début du gène (upstream)
   - Calcule la distance entre la fin du CDS et la fin du gène (downstream)
   - Gère le strand (+/-)

**Outputs** :
- `full_gff` : Annotations concaténées
- `full_fna` : Génomes concaténés
- `index` : Index FASTA
- `focal_proteins` : Protéome focal
- `local_db` : Base de données BLAST locale
- `geneCoords` : Coordonnées des gènes

---

### Étape 2 : Création des Maps Taxonomiques (`create_taxon_maps`)

**But** : Créer des dictionnaires pour le filtrage taxonomique

**Utilise** : `taxonkit` (outil de manipulation taxonomique NCBI)

**Outputs** :
- `strain2species.csv` : Mapping souche → espèce
- `eukaryotes.csv` : Liste des taxids eucaryotes (taxid 2759)

---

### Étape 3 : Recherche contre NR (`A2_NR.nf`)

**But** : Identifier les candidats potentiels avec des extensions manquantes

**Processus** :

1. **big_diamond** : Recherche Diamond contre NR (database eucaryotes taxid 2759)
   ```bash
   diamond blastp --query ${query} --db ${db} --taxonlist 2759 \
     --outfmt 6 qseqid sseqid qlen qstart qend qcovhsp scovhsp ppos staxids \
     --max-target-seqs 0 --fast -e 0.00001
   ```
   - Filtres : qcov > 60%, scov > 60%, ppos > 70%

2. **split_ids** : Division des IDs en N groupes (parallélisation)

3. **split_df** : Division du dataframe par groupe d'IDs

4. **parseDiamond** : Analyse des résultats pour identifier les candidats
   - Compte les espèces uniques
   - Calcule le % d'espèces avec alignement N-ter complet (qstart < 20)
   - Calcule le % d'espèces avec alignement C-ter complet (qlend < 20)
   
   **Critères de sélection** :
   - **N-ter candidats** : `percent_nter_q_15 <= 5%`
   - **C-ter candidats** : `percent_cter_q_15 <= 5%`

**Outputs** :
- `nter_candidates` : Liste des IDs candidats N-terminaux
- `cter_candidates` : Liste des IDs candidats C-terminaux
- `statistics` : Statistiques détaillées

---

### Étape 4 : Recherche d'Extensions (`A3_search_extensions.nf`)

**But** : Analyser en détail chaque candidat contre les espèces voisines locales

#### 4.1 Traitement N-terminal

**Flux** :
1. `grep_nter` : Extraction des protéines candidates
2. `ssearch_nter` : Alignement protéique contre le protéome local (ssearch36)
3. `parseNter` : Parsing et catégorisation

**Catégorisation** :
- **small_sstart** : Homologues avec sstart < 5 (proches du début)
  - Extension potentielle dans la région upstream du gène
- **big_sstart** : Homologues avec sstart >= 5 (éloignés du début)
  - L'alignement ne couvre pas le début de l'homologue

#### 4.2 Traitement C-terminal

**Flux similaire** :
- **small_send** : Homologues avec slend < 5
- **big_send** : Homologues avec slend >= 5

---

### Étape 5 : Analyse Nucléotidique (modules nter/ et cter/)

#### Pour les "small sstart/send" :

**Processus `small_elongate_and_align`** (script `bin/nter/small_nuc_search.sh`) :

1. **Élongation des sujets** :
   - Calcul de l'élongation : `elongation = qstart * 3 * 1.5` (arrondi multiple de 3)
   - Extension de la séquence génomique du sujet en amont (nter) ou aval (cter)
   - Utilisation de `bedtools slop` pour étendre les coordonnées

2. **Alignements multiples** :
   - `ssearch36` nucléotide vs nucléotide (séquence élongée)
   - `ssearch36` nucléotide vs nucléotide (séquence standard)
   - `ssearch36` protéine vs protéine (séquence élongée traduite)

3. **Parsing des alignements** (`parse_small_table`) :
   - Lecture du format "mA" de ssearch pour récupérer les alignements position par position
   - Comptage des gaps
   - Vérification des codons start/stop in-frame
   - Calcul de l'élongation récupérée

**Outputs par sujet** :
- `*_complete_small_ssearch_elong.tsv` : Résultats ssearch élongé
- `*_complete_small_ssearch_elong.aln` : Alignements détaillés
- `*_complete_small_ssearch_short.tsv` : Résultats ssearch standard
- `*_complete_small_tfastx.tsv` : Résultats protéiques
- `*_subjects_thresholds.tsv` : Seuils calculés par sujet

#### Pour les "big sstart/send" :

**Processus `ssearch_big_sstart/send`** (script `bin/nter/big_nuc_search.sh`) :

1. Extraction des CDS du sujet (sans élongation)
2. Alignement nucléotidique ssearch36
3. Comparaison de la position d'alignement avec la position protéique

---

## Modules Nextflow

### Modules Prétraitement (modules/preprocess/)

| Module | Fonction | Input | Output |
|--------|----------|-------|--------|
| `extract_protein` | Extraire protéome de GFF+FNA | (species, fna, gff) | *.faa |
| `makedb` | Créer base BLAST/Diamond | fasta | db directory / *.dmnd |
| `faidx` | Indexer génome | genome | *.fai |
| `concat` | Concaténer fichiers | files[] | concatenated.* |
| `sample_proteins` | Échantillonner protéines | (species, fna, gff), n | *.faa |
| `getGeneCoords` | Calculer coordonnées | gff | gene_coords.tsv |
| `checkAnnot` | Vérifier annotations uniques | - | - |

### Modules NR (modules/NR/)

| Module | Fonction | Input | Output |
|--------|----------|-------|--------|
| `big_diamond` | Recherche Diamond massive | query, db, taxid | *.tsv |
| `parseDiamond` | Parser résultats Diamond | df, strain2species, eukaryotes | nter_candidates, cter_candidates, stats |
| `split_ids` | Diviser IDs en N groupes | dataframe, n | *.txt files |
| `split_df` | Filtrer dataframe par IDs | df, ids | *.tsv |

### Modules d'Alignement (modules/)

| Module | Fonction | Input | Output |
|--------|----------|-------|--------|
| `blastp_strict` | BLASTp strict | query, db_dir | *.tsv |
| `blastp_souple` | BLASTp souple (short) | query, db_dir | *.tsv |
| `ssearch` | Alignement ssearch36 | query, subjects, cpus, mem | *.tsv |
| `seqkitGrep` | Extraire séquences par IDs | ids, fasta | output.fa |

### Modules N-terminal (modules/nter/)

| Module | Fonction |
|--------|----------|
| `ssearch_big_sstart` | Alignement pour gros sstart |
| `parse_big_table` | Parsing résultats gros sstart |
| `small_elongate_and_align` | Élongation + alignements multiples |
| `parse_small_table` | Parsing complexe avec analyse des gaps/codons |

### Modules C-terminal (modules/cter/)

| Module | Fonction |
|--------|----------|
| `ssearch_big_send` | Alignement pour gros send |
| `parse_big_table` | Parsing résultats gros send |
| `small_elongate_and_align` | Élongation + alignements multiples |
| `parse_small_table` | Parsing avec analyse des gaps/codons |

### Module Parsing Principal (`parseNeighborsBlastp.nf`)

**Processus `nter`** et **`cter`** : 
- Parsing Python avec Polars
- Filtrage homologie : evalue <= 1e-5, qcovhsp >= 70%, ppos > 50%
- Identification des "bad species" (espèces avec alignement complet)
- Catégorisation en big/small basée sur sstart/send
- Output en format Parquet par query

---

## Scripts Shell (bin/)

### bin/ssearch.sh
Wrapper pour ssearch36 avec parsing Polars :
```bash
ssearch36 -3 -p -s BL50 -f -11 -g -1 -T${ncpus} -XM${mem}G -m8BCL query subjects
```
- Matrice BL50, pénalités gap -11/-1
- Output format BLAST-like

### bin/nter/small_nuc_search.sh
Script principal pour l'analyse N-terminale des petits sstart :

1. **Lecture des paramètres** depuis fichier Parquet (via duckdb)
2. **Pour chaque sujet** :
   - Calcul élongation = max(qstart*3, upstream) * 1.5
   - Extraction CDS avec gffread
   - Extension avec bedtools slop (selon strand)
   - Extraction séquence élongée
   - Alignements multiples (ssearch nuc, ssearch prot, tfastx)
3. **Création des seuils** : threshold = elongation - 2

### bin/nter/big_nuc_search.sh
Script pour l'analyse N-terminale des gros sstart :
1. Conversion Parquet → TSV
2. Extraction CDS standards (sans élongation)
3. Alignement ssearch nucléotidique
4. Mapping qstart protéique → nucléotidique

### bin/cter/small_nuc_search.sh et big_nuc_search.sh
Équivalents pour l'analyse C-terminale (inversé : extension côté 3')

---

## Outils Utilisés

Tous les chemins d'outils sont **configurables** via `nextflow.config`. Par défaut, ils utilisent les noms d'outils qui doivent être dans le PATH.

### Outils Bioinformatiques

| Outil | Paramètre Config | Fonction |
|-------|------------------|----------|
| **Diamond** | `params.diamond` | Recherche similitude vs NR |
| **BLAST+** | `params.blastp`, `params.makeblastdb` | Alignements protéiques locaux |
| **FASTA36 (ssearch36)** | `params.ssearch` | Alignements Smith-Waterman |
| **tfasty36** | `params.tfasty` | Alignement prot vs ADN traduit |
| **gffread** | `params.gffread` | Extraction séquences depuis GFF |
| **bedtools** | `params.bedtools` | Manipulation coordonnées génomiques |
| **samtools** | `params.samtools` | Indexation génomes |
| **seqkit** | `params.seqkit` | Manipulation FASTA |
| **taxonkit** | `params.taxonkit` | Manipulation taxonomie NCBI |
| **AGAT** | conteneur `agat.sif` | Standardisation GFF |
| **faTrans** | `params.faTrans` | Traduction ADN → protéine |
| **faSize** | `params.faSize` | Taille séquences FASTA |
| **duckdb** | `params.duckdb` | Conversion Parquet → TSV |

### Librairies Python

```python
polars          # Traitement dataframes haute performance
gff3_parser     # Parsing fichiers GFF3
Bio (Biopython) # Manipulation séquences
csv, os, re, glob  # Utilitaires standard
```

---

## Configuration et Exécution

### Configuration (`nextflow.config`)

Le fichier de configuration centralise **tous les paramètres** et **chemins d'outils**. Les chemins peuvent être surchargés en ligne de commande ou via un fichier de configuration personnalisé.

```groovy
params {
    // Pipeline parameters
    n = 30                          // Nombre de groupes pour parallélisation
    output = "results/"             // Répertoire output
    
    // Input directories
    focal_dir = "input/focal/"      // Répertoire espèce focale
    neighbors_dir = "input/neighbors/"  // Répertoire espèces voisines
    
    // External databases
    nr = "/path/to/nr.dmnd"         // Base Diamond NR
    taxdump = "/path/to/taxdump.tar.gz"  // Données taxonomiques NCBI
    tmpdir = "/tmp"                 // Répertoire temporaire pour Diamond
    
    // --- OUTILS (tous configurables) ---
    // Alignment tools (FASTA36 suite)
    ssearch = "ssearch36"           // ou chemin absolu
    tfasty = "tfasty36"
    
    // Sequence manipulation
    seqkit = "seqkit"
    gffread = "gffread"
    faTrans = "faTrans"
    faSize = "faSize"
    
    // Database tools
    diamond = "diamond"
    blastp = "blastp"
    makeblastdb = "makeblastdb"
    
    // Taxonomy
    taxonkit = "taxonkit"
    
    // Data processing
    duckdb = "duckdb"
    
    // System tools
    samtools = "samtools"
    bedtools = "bedtools"
}

// Les outils sont exportés comme variables d'environnement pour les scripts shell
process {
    beforeScript = '''
        export SSEARCH="${params.ssearch}"
        export GFFREAD="${params.gffread}"
        # ... etc
    '''
}
```

### Surcharge des paramètres

```bash
# Via ligne de commande
nextflow run main.nf --ssearch /path/to/ssearch36 --diamond /path/to/diamond

# Via fichier de configuration personnalisé
nextflow run main.nf -c my_local.config
```

### Exécution

**Soumission PBS** :
```bash
qsub run.sh
# ou
qsub nxtflw_run.sh
```

**Exécution directe** :
```bash
nextflow run main.nf -resume
```

### Ressources HPC Typiques

| Queue | CPUs | RAM | Walltime | Utilisé pour |
|-------|------|-----|----------|--------------|
| bim | 70 | 400GB | 40000h | Diamond vs NR |
| bim | 16 | 120GB | 24h | Parsing BLAST |
| bim | 10 | 10GB | 15h | ssearch |
| common | 4-8 | 8-32GB | 24h | Prétraitement |
| lowprio | 4-6 | 4-8GB | variable | Parsing léger |

---

## Outputs

### Répertoire output/nter/
- `bigFinal` : Résultats consolidés des gros sstart
- `smallFinal` : Résultats consolidés des petits sstart
- `full.blast` : Tous les alignements filtrés
- `mRNA_species.tsv` : Mapping mRNA → espèce

### Répertoire output/cter/
- `bigFinal` : Résultats consolidés des gros send
- `smallFinal` : Résultats consolidés des petits send
- `full.blast` : Tous les alignements filtrés

### Colonnes Output (smallFinal nter)

| Colonne | Description |
|---------|-------------|
| qseqid | ID de la protéine query |
| sseqid | ID de l'homologue subject |
| gaps_query | Gaps dans la query |
| gaps_subject | Gaps dans le subject |
| qstart_nuc_elong | Position début alignement nucléotidique |
| sstart_nuc_gt_elongation | Si sstart > élongation attendue |
| stop_inframe | Codon STOP in-frame trouvé |
| start_inframe | Codon START in-frame trouvé |
| atg_on_elongated_subject_facing_query_start | ATG face au début query |
| meth_on_query_facing_subject_start | Méthionine query face au début subject |
| recovered_elongation | Ratio élongation récupérée |
| raw_recovered_elongation | Élongation en nucléotides |
| category | "small" ou "big" |

---

## Reconstruction du Pipeline

### Prérequis

1. **Installation Nextflow** :
```bash
curl -s https://get.nextflow.io | bash
mv nextflow ~/.local/bin/
```

2. **Installation des outils** (voir tableau ci-dessus)

3. **Configuration HPC PBS Pro**

### Étapes de Reconstruction

#### 1. Structure de base
```bash
mkdir -p DMEL/{flows,subworkflows/specifics,modules/{preprocess,NR,nter,cter,dev},bin/{nter,cter},src,container,input/{focal,neighbors},output/{nter,cter},test}
```

#### 2. Fichiers de configuration
- Créer `nextflow.config` avec les paramètres adaptés
- Créer `main.nf` qui appelle le workflow `align`

#### 3. Workflow Principal (flows/align.nf)
```groovy
include { preprocess_input_data } from '../subworkflows/A1_preprocess'
include { create_taxon_maps } from '../modules/create_taxon_maps'
include { NR } from '../subworkflows/A2_NR'
include { search_homologs } from '../subworkflows/A3_search_extensions'

workflow align {
    (full_gff, full_fna, index, focal_proteome, local_db, geneCoords, local_proteome) = preprocess_input_data()
    (strain2species, eukaryotes) = create_taxon_maps()
    (nter_candidates_IDs, cter_candidates_IDs) = NR(focal_proteome, strain2species, eukaryotes)
    search_homologs(nter_candidates_IDs, cter_candidates_IDs, local_db, focal_proteome, full_gff, full_fna, index, geneCoords, local_proteome)
}
```

#### 4. Subworkflows à implémenter
1. **A1_preprocess** : Enchaînement extraction → concat → makedb → faidx → getGeneCoords
2. **A2_NR** : Diamond → split → parseDiamond → collectFile
3. **A3_search_extensions** : grep → ssearch → parse → (big/small workflows)

#### 5. Modules clés à recréer

**Priorité 1** (flux principal) :
- `modules/preprocess/*.nf`
- `modules/NR/*.nf`
- `modules/blastp.nf`
- `modules/grepSeq.nf`
- `modules/parseNeighborsBlastp.nf`

**Priorité 2** (analyse détaillée) :
- `modules/nter/treatSmallSstart.nf` (le plus complexe)
- `modules/nter/treatBigSstart.nf`
- `modules/cter/treatSmallSend.nf`
- `modules/cter/treatBigSend.nf`

#### 6. Scripts shell critiques
- `bin/ssearch.sh` : Wrapper ssearch avec output formaté
- `bin/nter/small_nuc_search.sh` : Logique d'élongation et alignements multiples
- `bin/nter/big_nuc_search.sh` : Alignement sans élongation
- `bin/cter/*.sh` : Équivalents pour C-terminal

### Points d'attention

1. **Portabilité** : Tous les chemins d'outils sont désormais configurables dans `nextflow.config`. Les scripts shell utilisent des variables d'environnement avec fallback vers le PATH :
   ```bash
   SSEARCH=${SSEARCH:-ssearch36}  # Utilise $SSEARCH si défini, sinon 'ssearch36'
   ```

2. **Format Parquet** : Les données intermédiaires utilisent Parquet, nécessitant duckdb pour la conversion

3. **Parsing ssearch** : Le format `-m A` de ssearch36 produit un alignement position par position qui nécessite un parsing spécifique

4. **Gestion des strands** : L'élongation doit être faite du bon côté selon le strand (+/-)

5. **Seuils et filtres** :
   - Homologie : evalue <= 1e-5, qcovhsp >= 70%, ppos > 50%
   - Candidats NR : < 5% avec alignement complet aux extrémités
   - Catégorisation : sstart/send < 5 vs >= 5

---

## Licence et Contact

*Pipeline développé par Simon Herman - BIM Team*

---

## Changelog

- **v1.0** : Pipeline initial pour *D. melanogaster*
- Support pour analyse N-terminale et C-terminale
- Intégration Diamond, BLAST, ssearch36
- Parallélisation PBS Pro

