#!/usr/bin/bash
# big_nuc_search.sh - Big sstart nucleotide search for N-terminal extension analysis
# Tools paths can be set via environment variables or will use PATH defaults

query_and_subjects=$1
fna=$2
gff=$3

# Use environment variables if set, otherwise use tool names (expects them in PATH)
SSEARCH=${SSEARCH:-ssearch36}
GFFREAD=${GFFREAD:-gffread}
DUCKDB=${DUCKDB:-duckdb}

echo "" > log.txt

is_integer() {

     [[ $1 =~ ^[0-9]+$ ]]

}


if [ -s "$query_and_subjects" ]; then

    # For some reasons, the previous process is not able to correctly write the file in plain text
    # We wrote it in parquet ( binary ) then convert it to tsv here

    $DUCKDB -c "COPY (SELECT * FROM read_parquet('$query_and_subjects')) TO 'query_subjects.tsv' (FORMAT 'CSV', HEADER 'False', DELIMITER '\t');"

    # extract query before entering the loop. Assume that there is no header
    query=$(awk -F "\t" 'NR==1 {print $1}' query_subjects.tsv)
    query_pattern=$(echo "${query}" | sed 's/\./\\./g' | sed 's/|/\\|/g')

    echo "Query : " >> log.txt
    echo "$query" >> log.txt
    echo "$query_pattern" >> log.txt

    grep -E "${query_pattern}\b" "${gff}" > query_cds
    cut -f2 query_subjects.tsv | sort -u > subject_ids.txt
    grep -F -f subject_ids.txt "${gff}" > subjects_gff_subset
    $GFFREAD -x - -g "${fna}" query_cds > query.fna

    echo -e "sseqid\tqstart\tqstart_mapped" > "${query}_qstarts.tsv"

    while IFS= read -r line; do
        
        IFS=$'\t' read -r -a fields <<< "$line"

        # Structure of line : query, subject, qstart, sstart, upstream
        # tab separated 

        subject=${fields[1]}

        # Qstart and sstart are from the proteic alignment output, so they are relative to a protein sequence
        qstart=${fields[2]}
        sstart=${fields[3]}

        qstart_nuc_mapped=$((qstart * 3))

        if ! is_integer "$qstart" || ! is_integer "$sstart" ; then
            echo "One of the values is not an integer : " >> log.txt
            echo "$qstart" >> log.txt
            echo "$sstart" >> log.txt
            echo "Exiting" >> log.txt
            exit 1 
        fi

        echo "Subject : " >> log.txt
        echo "$subject" >> log.txt

        subject_pattern=$(echo "$subject" | sed 's/\./\\./g' | sed 's/|/\\|/g')

        echo "Subject pattern : " >> log.txt
        echo "$subject_pattern" >> log.txt

        rm -f subject_cdss

        subject_pattern=$(echo "$subject" | sed 's/\./\\./g' | sed 's/|/\\|/g')
        grep -E "$subject_pattern\b" subjects_gff_subset | $GFFREAD -C | grep "CDS" | sort -n -k4 > subject_cdss

        if [ -s subject_cdss ]; then

            $GFFREAD -x - -g "${fna}" subject_cdss >> subjects.fna

        else 

            echo "No CDS found for $subject" >> log.txt
            exit 1

        fi

        echo -e "$subject\t$qstart\t$qstart_nuc_mapped" >> "${query}_qstarts.tsv"

    done < "query_subjects.tsv"

    echo -e "qseqid\tqlen_nuc\tsseqid\tslen_nuc\tpident\talign_length\tbar\tgaps\tqstart_nuc\tqend_nuc\tsstart_nuc\tsend_nuc\tevalue\tbitscore\tCIGAR" > "${query}_big_ssearch.tsv"

    $SSEARCH -3 -n -T8 -XM8G -E 10000 -m8BCL query.fna subjects.fna >> "${query}_big_ssearch.tsv"

else

    touch "${query}_big_ssearch.tsv"
    touch "${query}_subjects_thresholds.tsv"
    echo "No query_subjects.tsv found" >> log.txt

fi