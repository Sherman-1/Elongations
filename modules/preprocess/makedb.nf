process blast_makedb {


    clusterOptions '--partition=bim --cpus-per-task=4 --mem=8gb'

    input:
    path fasta

    output:
    stdout

    script:
    """
    #!/usr/bin/bash

    ${params.makeblastdb} -in $fasta -dbtype 'prot' -out blast -parse_seqids -hash_index > /dev/null
    echo "\$PWD"

    """

    
}

process diamond_makedb {


    clusterOptions '--partition=bim --cpus-per-task=8 --mem=30gb --time=1000:00:00 --nodelist=node04'

    input:
    path fasta

    output:
    path "${fasta.baseName}"

    script:
    """
    #!/usr/bin/bash

    ${params.diamond} makedb --in $fasta --db "${fasta.baseName}.dmnd" --threads 8
    """

    
}