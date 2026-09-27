process small_diamond {

    clusterOptions '--partition=bim --cpus-per-task=30 --mem=80gb --time=1000:00:00'

    label 'small_diamond'

    input:
    path query
    path db
    val taxid

    output:
    path "*.tsv"

    script:
    """
    #!/usr/bin/bash

    output=\$(basename "${query}" | cut -d. -f1)

    DIAMOND=${params.diamond}

    if [ ${taxid} -gt 0 ]; then

        \$DIAMOND blastp --query ${query} --db ${db} \\
            --outfmt 6 qseqid sseqid qlen slen qstart qend sstart send evalue qcovhsp ppos scovhsp staxids \\
            --threads 30 --taxonlist ${taxid} --header simple \\
            -k0 --verbose --sensitive \\
            -b 5 -c 1 --out "${query.baseName}.tsv"

    else
        \$DIAMOND blastp --query ${query} --db ${db} \\
            --outfmt 6 qseqid sseqid qlen slen qstart qend sstart send evalue qcovhsp ppos scovhsp \\
            --threads 30 --header simple \\
            -k0 --verbose --sensitive \\
            -b 5 -c 1 --out "${query.baseName}.tsv"

    fi
    """
}

process big_diamond {

    clusterOptions '--partition=bim --cpus-per-task=70 --mem=400gb --time=248:00:00 --nodelist=node04'

    label 'big_diamond'

    input:
    path query
    path db
    val taxid

    output:
    path "*.tsv"

    script:
    """
    #!/usr/bin/bash

    DIAMOND=${params.diamond}
    TMPDIR=${params.tmpdir}

    echo "Running DIAMOND on ${query} with taxid ${taxid}"

    if [ ${taxid} -gt 0 ]; then

        echo -e "qseqid\\tsseqid\\tqlen\\tqstart\\tqlend\\tstaxids" > "${query.baseName}.tsv"

        \$DIAMOND blastp --query ${query} --db ${db} --tmpdir \$TMPDIR \\
            --outfmt 6 qseqid sseqid qlen qstart qend qcovhsp scovhsp ppos staxids \\
            --threads 70 --taxonlist ${taxid} \\
            -k0 --verbose \\
            --sensitive -b 15 -c 1 -e 0.00001 --out temp

        awk 'BEGIN {OFS = FS = "\\t";} \$6 > 60 && \$7 > 60 && \$8 > 70 { print \$1,\$2,\$3,\$4,\$3-\$5,\$9}' temp >> "${query.baseName}.tsv"
        # Filter on qcov, scov and ppos; only keep the columns we need later
        rm temp

    else

        \$DIAMOND blastp --query ${query} --db ${db} --tmpdir \$TMPDIR \\
            --outfmt 6 qseqid sseqid qlen slen qstart qend sstart send evalue qcovhsp scovhsp ppos \\
            --threads 70 --header simple \\
            -k0 --verbose \\
            --sensitive -b 15 -c 1 -e 0.00001 --out "${query.baseName}.tsv"

    fi
    """
}
