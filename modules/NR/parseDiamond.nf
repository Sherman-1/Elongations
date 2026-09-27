process parseDiamond {

    clusterOptions '--partition=bim --cpus-per-task=6 --mem=50gb --time=1000:00:00'    


    input:
    each path(df)
    path strain2species
    path eukaryotes

    output:
    path 'nter_candidates_*.tsv', emit : nter
    path 'cter_candidates_*.tsv', emit : cter
    path 'statistics_*.tsv', emit : stats


    script:
    """
    #!/usr/bin/bash
    python3 ${projectDir}/bin/parse_diamond.py ${df} ${strain2species} ${eukaryotes}
    """
    
}
