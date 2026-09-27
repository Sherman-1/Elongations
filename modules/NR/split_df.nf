process split_df {

    clusterOptions '--partition=kill-requeue --cpus-per-task=6 --mem=8gb --time=48:00:00'


    input:
    path df
    each path(ids)

    output:
    path "*.tsv", emit : divided

    script:
    """
    #!/usr/bin/env bash

    split_dataframe.sh $df $ids
    """

}
