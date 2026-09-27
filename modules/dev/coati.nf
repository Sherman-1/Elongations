process nter {

    clusterOptions '--partition=lowprio --cpus-per-task=1 --mem=2gb'

    errorStrategy { sleep(Math.pow(2, task.attempt) as long); return 'retry' }
    maxRetries 10
    
    label 'nter'

    input:
    path candidate_and_homolog
    path gff
    path fna 

    output:
    path "*.result"
    path "*.tsv"

    script:
    """
    #!/usr/bin/bash


    coati.sh ${candidate_and_homolog} ${gff} ${fna} "forward"

    """

}

process cter {

    clusterOptions '--partition=lowprio --cpus-per-task=1 --mem=2gb --time=24:00:00'

    errorStrategy { sleep(Math.pow(2, task.attempt) * 1 as long); return 'retry' }
    maxRetries 20
    
    label 'cter'

    input:
    path candidate_and_homolog
    path gff
    path fna 

    output:
    path "*.tsv"

    script:
    """
    #!/usr/bin/bash


    coati.sh ${candidate_and_homolog} ${gff} ${fna} "reverse"

    """

}

