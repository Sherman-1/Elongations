process sample_proteins { 

    label 'sample'

    clusterOptions '--partition=common --cpus-per-task=4 --mem=6gb --time=24:00:00'

    input:
    tuple val(species), path(fna), path(gff)
    val number

    output:
    path "*.faa"
    script:
    """
    #!/usr/bin/bash
    
    grep -F -f <(gffread -C ${gff} | grep "CDS" | cut -f 9,9 | sort -u | shuf -n ${number}) ${gff} \
    | gffread -J -y - -g ${fna} /dev/stdin | seqkit seq --remove-gaps --max-len 1500 - > ${species}.faa
    
    """
}