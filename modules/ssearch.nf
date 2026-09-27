process ssearch {

    clusterOptions '--partition=kill-requeue --cpus-per-task=8 --mem=8gb --time=24:00:00'

    label 'ssearch'

    input:
    path candidate_and_homologs
    path gff
    path fna 
    path index
    
    output:
    path "*_ssearch.tsv"
    path "*_ssearch_aligns.txt"

    script:
    """
    
    ssearch.sh ${candidate_and_homologs} ${fna} ${gff} ${index}

    """

}