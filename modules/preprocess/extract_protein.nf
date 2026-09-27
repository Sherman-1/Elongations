process extract_protein {

    clusterOptions '--partition=common --cpus-per-task=4 --mem=8gb --time=24:00:00'

    label "extract_protein"

    input:
    tuple val(species), path(fna), path(gff)

    output:
    path "*.faa"

    script:
    """
    #!/usr/bin/bash

    gffread -J -y - -g ${fna} ${gff} | seqkit seq --remove-gaps - > ${species}.faa

    """
}
