process small_elongate_and_align { 

    clusterOptions '--partition=bim --cpus-per-task=10 --mem=10gb --time=10:00:00'

    input:
    each path(query_subjects)
    path genomes
    path annotations
    path index

    output:
    tuple path("*_complete_small_ssearch_elong.tsv"), 
          path("*_complete_small_ssearch_elong.aln"), 
          path("*_subjects_thresholds.tsv"), 
          path("elongated_subjects.fna"), 
          path("elongated_subjects.faa"),
          path("*_complete_small_tfastx.tsv"), 
          path("*_complete_small_ssearch_short.tsv")

    script:
    """
    #!/usr/bin/bash
    ${projectDir}/bin/cter/small_nuc_search.sh ${query_subjects} ${genomes} ${annotations} ${index}
    """

}


process parse_small_table {

    clusterOptions '--partition=kill-requeue --cpus-per-task=4 --mem=4gb --time=2:00:00'

    input:
    tuple path(complete_small_ssearch_elong_table), 
          path(complete_small_ssearch_elong_aln), 
          path(thresholds), 
          path(elongated_subjects_fna), 
          path(elongated_subjects_faa),
          path(complete_small_tfastx_table),
          path(complete_small_ssearch_short_table)

    output:
    path "*_final_res.tsv", emit : small_send_res

    script:
    """
    #!/usr/bin/bash
    # rerun 2026-08-08 : -E 10000
    python3 ${projectDir}/bin/cter/parse_small_table.py \
        ${complete_small_ssearch_elong_table} \
        ${complete_small_ssearch_elong_aln} \
        ${complete_small_ssearch_short_table} \
        ${complete_small_tfastx_table} \
        ${thresholds} \
        ${elongated_subjects_faa}
    """

}
