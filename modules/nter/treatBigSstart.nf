process ssearch_big_sstart { 

    clusterOptions '--partition=bim --cpus-per-task=10 --mem=10gb --time=15:00:00'

    input:
    each path(query_subjects)
    path genomes
    path annotations

    output:
    tuple path("*_big_ssearch.tsv"),
          path("*_qstarts.tsv")
    
    script:
    """
    #!/usr/bin/bash

    ${projectDir}/bin/nter/big_nuc_search.sh ${query_subjects} ${genomes} ${annotations}
    """
}


process parse_big_table { 

    clusterOptions '--partition=bim --cpus-per-task=4 --mem=4gb --time=15:00:00'

    input:
    tuple path(ssearch_big_table),
          path(qstarts)


    output:
    path "*_final_res.tsv", emit : big_sstart_res

    script:
    """
    #!/usr/bin/bash
    python3 ${projectDir}/bin/nter/parse_big_table.py ${ssearch_big_table} ${qstarts}
    """

}