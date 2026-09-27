
process nter {


    clusterOptions '--partition=bim --cpus-per-task=16 --mem=120gb --time=24:00:00'

    publishDir "output/nter/", mode : 'copy', pattern : "*.blast", overwrite : true
    publishDir "output/nter/", mode : 'copy', pattern : "mRNA_species.tsv", overwrite : true
    publishDir "output/nter/", mode : 'copy', pattern : "bad_species.tsv", overwrite : true

    label 'nter'

    input:
    path candidates_vs_neighbors_blastp_output
    path coding_coords

    output:
    path "*_small.parquet", emit : small 
    path "*_big.parquet", emit : big 
    path "bad_species.tsv", emit : bad_species
    path "full.blast", emit : full_blast



    script:
    """
    #!/usr/bin/bash

    echo "Parsing blastp output for N-terminal candidates"
    echo "Ca repart les amis!!"
    python3 ${projectDir}/bin/parse_blastp_nter.py ${candidates_vs_neighbors_blastp_output} ${coding_coords} ${projectDir}/input/neighbors/
    """

}

process cter {


    clusterOptions '--partition=bim --cpus-per-task=16 --mem=120gb --time=24:00:00'

    publishDir "output/cter/", mode : 'copy', pattern : "*.blast", overwrite : true
    publishDir "output/cter/", mode : 'copy', pattern : "mRNA_species.tsv", overwrite : true
    publishDir "output/cter/", mode : 'copy', pattern : "bad_species.tsv", overwrite : true

    label 'cter'

    input:
    path candidates_vs_neighbors_blastp_output
    path coding_coords

    output:
    path "*_small.parquet", emit : small 
    path "*_big.parquet", emit : big 
    path "bad_species.tsv", emit : bad_species
    path "full.blast", emit : full_blast



    script:
    """
    #!/usr/bin/bash
    echo "Parsing blastp output for C-terminal candidates, new run with scovhsp >= 70"
    python3 ${projectDir}/bin/parse_blastp_cter.py ${candidates_vs_neighbors_blastp_output} ${coding_coords} ${projectDir}/input/neighbors/
    """

}
