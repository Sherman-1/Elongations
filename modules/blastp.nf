process blastp_strict {

    clusterOptions '--partition=bim --cpus-per-task=20 --mem=40gb --time=10:00:00'

    input:
    path query
    path db_dir

    output:
    path "*.tsv"

    script:
    """
    #!/bin/bash

    echo -e "qseqid\\tsseqid\\tqlen\\tqstart\\tqend\\tsstart\\tsend\\tslen\\tevalue\\tqcovhsp\\tppos" > "${query.baseName}_vs_neighbors.tsv"
    ${params.blastp} -query ${query} -db ${db_dir}/blast -num_threads 20 -max_hsps 1 -evalue 0.0001 \\
    -outfmt '6 qseqid sseqid qlen qstart qend sstart send slen evalue qcovhsp ppos' >> "${query.baseName}_vs_neighbors.tsv"

    """


}

process ssearch { 

    input:
    path query
    path subjects
    val number_cpus
    val mem 

    clusterOptions { "--partition=bim --time=480:00:00 --mem=${mem} --cpus-per-task=${number_cpus}" }

    output:
    path "*.tsv"

    script:
    """
    #!/bin/bash
    # Proteic alignment, BP50 matrix, -11/-1 gap penalties 
    ${projectDir}/bin/ssearch.sh ${query} ${subjects} ${number_cpus} ${mem}   
    """
}
