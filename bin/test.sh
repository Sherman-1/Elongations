#!/bin/bash
# test.sh - Test script for debugging
# Tools paths can be set via environment variables or will use PATH defaults

query_and_subjects=$1
fna=$2
gff=$3
index=$4

# Use environment variables if set, otherwise use tool names (expects them in PATH)
SSEARCH=${SSEARCH:-ssearch36}
GFFREAD=${GFFREAD:-gffread}

echo "" > log.txt
cut -f1,2 "$index" > neighbors.genome

# Candidates and homologs is a tsv, query is the first column and is always the same
# extract query before entering the loop

query=$(awk 'BEGIN { FS = "\t";}' 'FNR == 1 { print $1 }' "${query_and_subjects}")
query_pattern=$(echo "${query}" | sed 's/\./\\./g' | sed 's/|/\\|/g')

echo "$query" >> log.txt
echo "$query_pattern" >> log.txt

grep -E "${query_pattern}\b" "${gff}" > query_cds
$GFFREAD -x - -g "${fna}" query_cds > query.fna

echo "query fasta : " >> log.txt
cat query.fna >> log.txt


while IFS= read -r line; do
    
    IFS=$'\t' read -r -a fields <<< "$line"

    # Structure of line : query, subject, qstart, sstart, upstream
    # tab separated 

    subject=${fields[2]}
    upstream=${fields[3]}

    echo "Subject : " >> log.txt
    echo "$subject" >> log.txt

    subject_pattern=$(echo "$subject" | sed 's/\./\\./g' | sed 's/|/\\|/g')

   
    # Make sure that the elongation is a multiple of 3
    # and is based on the length of the elongation
    buffer=$(echo "scale=0; (($upstream * 1.5) + 0.5) / 1" | bc)
    remainder=$(echo "$buffer % 3" | bc)

    if [ "$remainder" -eq 0 ]; then 
        elongation="$buffer"
    else 
        elongation="$(echo "$buffer + (3 - $remainder)" | bc)"
    fi
    
    echo "upstream : $upstream" >> log.txt 
    echo "elongation : $elongation" >> log.txt

    rm -f subject_cdss
    rm -f modified_coords

    
    subject_pattern=$(echo "$subject" | sed 's/\./\\./g' | sed 's/|/\\|/g')
    grep -E "$subject_pattern\b" "${gff}" | $GFFREAD -C | grep "CDS" | sort -n -k4 > subject_cdss

    exons_number=$(wc -l subject_cdss | awk '{print $1}')

    if [ $exons_number -eq 1 ]; then

		bedtools slop -i subject_cdss -g neighbors.genome -b $elongation >> modified_coords
        
									     
    else

		rm -f modified_coords

		exon_counter=1

		while read -r exon; do

			if [ $exon_counter -eq 1 ]; then

				echo "$exon" | bedtools slop -i - -g neighbors.genome -l $elongation -r 0 >> modified_coords

			elif [ $exon_counter -eq $exons_number ]; then

				echo "$exon" | bedtools slop -i - -g neighbors.genome -l 0 -r $elongation >> modified_coords

			else

				echo "$exon" >> modified_coords

			fi

			exon_counter=$((exon_counter+1))

		done < subject_cdss

    fi

    if [ -s modified_coords ]; then

        $GFFREAD -x - -g "${fna}" modified_coords >> subjects.fna

    else

        echo "No cds found for $subject" >> log.txt
        exit 1

    fi
    

done < "$query_and_subjects"


# Perform the search
echo -e "qseqid\tsseqid\tppos\tqlen\tmismatch\tgaps\tqstart\tqend\tsstart\tsend\tevalue\tbitscore" > "${query}_big_ssearch.tsv"
$SSEARCH -3 -T 6 -n -m8 query.fna subjects.fna >> "${query}_big_ssearch.tsv"
$SSEARCH -3 -T 6 -n query.fna subjects.fna >> "${query}_big_ssearch_aligns.txt"
