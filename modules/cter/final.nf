process analyze_small_cter { 

    clusterOptions '--partition=common --cpus-per-task=16 --mem=32gb --time=15:00:00'

    input:
    path big_cter


    output:
    path "*_analyzed.tsv", emit : analyzed

    script:
    """
    #!/usr/bin/bash
    python3 ${projectDir}/bin/cter/analyze_small.py ${big_cter}
    """

}