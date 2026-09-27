process faidx {

    clusterOptions '--partition=common --cpus-per-task=4 --mem=6gb --time=1000:00:00'

    input:
    file genome 

    output:
    file "${genome}.fai"

    script:
    """
    samtools faidx ${genome} > ${genome}.fai
    """

}