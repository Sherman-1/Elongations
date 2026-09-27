process checkAnnot {

    clusterOptions '--partition=lowprio --cpus-per-task=4 --mem=64gb --time=10:00:00'


    script:
    """
    #!/usr/bin/bash
    python3 ${projectDir}/bin/check_annotation.py "input/*/*.gff"
    """

}