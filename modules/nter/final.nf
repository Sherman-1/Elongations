process analyze_small_nter { 

    clusterOptions '--partition=common --cpus-per-task=16 --mem=32gb --time=15:00:00'

    input:
    path small_nter


    output:
    path "*_analyzed.tsv", emit : analyzed

    script:
    """
    #!/usr/bin/bash
    python3 ${projectDir}/bin/nter/analyze_small.py ${small_nter}
    """

}