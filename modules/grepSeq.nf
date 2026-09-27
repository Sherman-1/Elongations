process seqkitGrep {


    clusterOptions '--partition=kill-requeue --cpus-per-task=4'

    input:
    path ids
    path fasta

    output:
    path "output.fa"    

    script:
    """
    ${params.seqkit} grep -f ${ids} ${fasta} > output.fa
    """


}

process blastdbcmd { 

    clusterOptions '--partition=kill-requeue --cpus-per-task=4'

    input:
    path db
    path ids

    output:
    path "output.fa"    

    script:
    """

    blastdbcmd -db ${db} -entry ${ids} -out output.fa
    """
}