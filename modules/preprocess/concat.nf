process concat { 

    clusterOptions '--partition=common --cpus-per-task=4 --mem=6gb --time=24:00:00'


    input:
    path files

    output:
    path "concatenated.*"

    script:
    """
    filetype=\$(basename ${files[0]} | rev | cut -d. -f1 | rev)

    cat $files > concatenated.\${filetype}
    """
}