process getGeneCoords {

    clusterOptions '--partition=bim --cpus-per-task=16 --mem=64gb --time=10:00:00'

    input:
    path gff

    output:
    path "*.tsv" , emit : geneCoords

    script:
    """
    #!/usr/bin/bash
    python3 ${projectDir}/bin/get_gene_coords.py ${gff}
    """

}
