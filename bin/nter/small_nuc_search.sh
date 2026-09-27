#!/bin/bash
# small_nuc_search.sh - Small sstart nucleotide search for N-terminal extension analysis
# Tools paths can be set via environment variables or will use PATH defaults

query_and_subjects=$1
fna=$2
gff=$3
index=$4

# Use environment variables if set, otherwise use tool names (expects them in PATH)
SSEARCH=${SSEARCH:-ssearch36}
TFASTY=${TFASTY:-tfasty36}
GFFREAD=${GFFREAD:-gffread}
FATRANS=${FATRANS:-faTrans}
DUCKDB=${DUCKDB:-duckdb}
FASIZE=${FASIZE:-faSize}
SEQKIT=${SEQKIT:-seqkit}
BEDTOOLS=${BEDTOOLS:-bedtools}

echo "" > log.txt

is_integer() {

     [[ $1 =~ ^[0-9]+$ ]]

}

extract_alignment_section() {

    local input_file="$1"
    local output_file="$2"


    if [[ ! -f "$input_file" ]]; then
        echo "Input file does not exist." >&2
        return 1
    fi

    rm -f temp_file.txt

    # Non empty lines
    awk '
    NF {print > "temp_file.txt"} 
    ' "$input_file"

    awk '
    BEGIN {start=0}
    /^>>/ {start=1}
    /^[0-9]+ residues in [0-9]+ query\s+sequences$/ {exit}
    /^[0-9]+ residues in [0-9]+ library\s+sequences$/ {exit}
    start {print}
    ' temp_file.txt >> "$output_file"

    rm temp_file.txt
}

echo "Bonjour" >> log.txt


if [ -s "$query_and_subjects" ]; then

    # For some reasons, the previous process is not able to correctly write the file in plain text
    # We write it in parquet ( binary ) then convert it to tsv here
    $DUCKDB -c "COPY (SELECT * FROM read_parquet('$query_and_subjects')) TO 'query_subjects.tsv' (FORMAT 'CSV', HEADER 'False', DELIMITER '\t');"

    cut -f1,2 "$index" > neighbors.genome

    # Candidates and homologs is a tsv representing proteic alignments between a query and its subjects
    # Columns of df : query, subject, qstart, sstart, upstream, strand
    # upstream is the distance between the start of the subject protein and the start its corresponding gene ( kind of UTR )


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

    : > "${query}_subjects_thresholds.tsv"
    : > elongated_subjects.fna
    : > standard_subjects.fna
    : > "${query}_fragment_small_ssearch_elong.aln"

    echo -e "qseqid\tqlen_nuc_elong\tsseqid\tslen_nuc_elong\tpident\talign_length\tbar\tgaps\tqstart_nuc_elong\tqend_nuc_elong\tsstart_nuc_elong\tsend_nuc_elong\tevalue\tbitscore\tCIGAR" > "${query}_fragment_small_ssearch_elong.tsv"
    echo -e "qseqid\tqlen_nuc_short\tsseqid\tslen_nuc_short\tpident\talign_length\tbar\tgaps\tqstart_nuc_short\tqend_nuc_short\tsstart_nuc_short\tsend_nuc_short\tevalue\tbitscore\tCIGAR" > "${query}_fragment_small_ssearch_short.tsv"
    echo -e "qseqid\tqlen_tfastx\tsseqid\tslen_tfastx\tpident\talign_length\tbar\tgaps\tqstart_tfastx\tqend_tfastx\tsstart_tfastx\tsend_tfastx\tevalue\tbitscore\tCIGAR" > "${query}_fragment_small_tfastx.tsv"



    while IFS= read -r line; do
        
        IFS=$'\t' read -r -a fields <<< "$line"

        # Structure of line : query, subject, qstart, sstart, upstream
        # tab separated 
        subject=${fields[1]}

        # Qstart and sstart are from the proteic alignment output, so they are relative to a protein sequence
        qstart=${fields[2]}
        sstart=${fields[3]}
        

        strand=${fields[5]}
        upstream=${fields[4]}

        # These two are not meant to point the exact nucleotide that gave the amino acid at position sstart or qstart
        # It is rather the portion of the query and subject that does not align in the Nter region, multiplied by 3
        # to give its length in nucleotides
        qstart_mapped_nuc=$(((qstart - 1) * 3))
        sstart_mapped_nuc=$(((sstart - 1) * 3))

        if ! is_integer "$qstart" || ! is_integer "$sstart" || ! is_integer "$upstream"; then
            echo "One of the values is not an integer : " | tee -a log.txt >&2
            echo "$qstart" >> log.txt
            echo "$sstart" >> log.txt
            echo "$upstream" >> log.txt
            echo "Exiting" >> log.txt
            exit 1 
        fi

        echo "Subject : " >> log.txt
        echo "$subject" >> log.txt

        subject_pattern=$(echo "$subject" | sed 's/\./\\./g' | sed 's/|/\\|/g')

        echo "Subject pattern : " >> log.txt
        echo "$subject_pattern" >> log.txt

        # We want the elongation to be the greatest between qstart_nuc and upstream. The goal here 
        # is to not miss any potential upstream region when screening via nucleotides alignments

        ### /!\ /!\ /!\ 

        # The code is modified for now, we only keep the length of qstart x2, sometimes 
        # upstream elongates the subject gene too far and matches are sometimes found 
        # in this upper region which parasitizes the results

        ### /!\ /!\

        if [ $qstart_mapped_nuc -gt $upstream ]; then

            temp=$qstart_mapped_nuc

        else
        
            temp=$qstart_mapped_nuc # To bring back the old way, just change this line to temp=$upstream
        fi

        echo "query_subject_thresholds at this point : " >> log.txt
        cat "${query}_subjects_thresholds.tsv" >> log.txt

        # We want the elongation to be 1.5 times the value of temp, and then we want to round it up to the nearest multiple of 3
        # Euclidean division by 3, if the remainder is 0, we keep the value, otherwise we add the remainder to the value
        buffer=$(echo "scale=0; (($temp * 1.5)) / 1" | bc)
        remainder=$(echo "$buffer % 3" | bc)
        if [ "$remainder" -eq 0 ]; then 
            elongation="$buffer"
        else 
            elongation="$(echo "$buffer + (3 - $remainder)" | bc)"
        fi

        echo "qstart_nuc : $qstart_mapped_nuc" >> log.txt
        echo "sstart_nuc : $sstart_mapped_nuc" >> log.txt
        echo "upstream : $upstream" >> log.txt 
        echo "elongation : $elongation" >> log.txt
        echo "strand : $strand" >> log.txt

        rm -f standard_subject_coords
        rm -f modified_subject_coords

        grep -E "$subject_pattern\b" subjects_gff_subset | $GFFREAD -C | grep "CDS" | sort -n -k4 > standard_subject_coords

        echo "Standard coords for $subject" >> log.txt
        cat standard_subject_coords >> log.txt

        exons_number=$(wc -l standard_subject_coords | awk '{print $1}')

        if [ $exons_number -eq 1 ]; then

            echo "Found only one exon for $subject" >> log.txt

            if [ "$strand" == "+" ]; then

                $BEDTOOLS slop -i standard_subject_coords -g neighbors.genome -l $elongation -r 0 > modified_subject_coords

            elif [ "$strand" == "-" ]; then

                $BEDTOOLS slop -i standard_subject_coords -g neighbors.genome -l 0 -r $elongation > modified_subject_coords

            fi
            
                                                
        else 

            echo "Found multiple exons for $subject" >> log.txt

            exon_counter=1

            while read -r exon; do
                
                # First exon of the mRNA
                if [ $exon_counter -eq 1 ]; then


                    # If the subject gene is on the positive strand, we elongate it " on the left "
                    # hence we add $elongation to the left part of the first exon
                    # If not, we do not touch the first exon ( as it is the end of the gene )
                    if [ "$strand" == "+" ]; then

                        echo "$exon" | $BEDTOOLS slop -i - -g neighbors.genome -l $elongation -r 0 >> modified_subject_coords

                    elif [ "$strand" == "-" ]; then

                        echo "$exon" >> modified_subject_coords

                    fi

                # Last exon of the mRNA
                elif [ $exon_counter -eq $exons_number ]; then

                    # If the subject gene is on the negative strand, we elongate it " on the right "
                    # hence we add $elongation to the right part of the last exon
                    # If not, we do not touch the last exon ( as it is the end of the gene )
                    if [ "$strand" == "-" ]; then

                        echo "$exon" | $BEDTOOLS slop -i - -g neighbors.genome -l 0 -r $elongation >> modified_subject_coords

                    elif [ "$strand" == "+" ]; then

                        echo "$exon" >> modified_subject_coords

                    fi

                # The exons in the middle are not touched
                else

                    echo "$exon" >> modified_subject_coords

                fi

                exon_counter=$((exon_counter+1))

                

            done < standard_subject_coords

        fi

        if [ -s modified_subject_coords ]; then

            echo "File modified_subject_coords is not empty" >> log.txt

            # For each subject, we subseq the query to only take
            # the n-ter segment that does not align with the subject ( + a little buffer )
            query_fragment=$(echo "$qstart_mapped_nuc * 1.5" | bc | cut -d'.' -f1)
            echo "Query fragment : $query_fragment" >> log.txt

            echo "Extracting query fragment" >> log.txt
            $SEQKIT subseq -r 1:$query_fragment query.fna > query_fragment.fna

            echo "Extracting subject $subject sequences" >> log.txt
            $GFFREAD -x - -g "${fna}" modified_subject_coords > temp_elongated_subject.fna
            $GFFREAD -x - -g "${fna}" standard_subject_coords > temp_standard_subject.fna

            $FATRANS query_fragment.fna query_fragment.faa
            $FATRANS temp_elongated_subject.fna temp_elongated_subject.faa

            $SSEARCH -3 -n -T8 -XM8G -E 10000 -m8BCL -m "FA fragment_align_elong.temp" query_fragment.fna temp_elongated_subject.fna >> "${query}_fragment_small_ssearch_elong.tsv"
            $SSEARCH -3 -n -T8 -XM8G -E 10000 -m8BCL query_fragment.fna temp_standard_subject.fna >> "${query}_fragment_small_ssearch_short.tsv"
            $SSEARCH -3 -p -T8 -XM8G -E 10000 -m8BCL -s BL50 -f -11 -g -1 query_fragment.faa temp_elongated_subject.faa >> "${query}_fragment_small_tfastx.tsv"

            cat temp_elongated_subject.fna >> elongated_subjects.fna
            cat temp_standard_subject.fna >> standard_subjects.fna

            extract_alignment_section fragment_align_elong.temp "${query}_fragment_small_ssearch_elong.aln"

            # Fasize -detailed : id \t length
            elongated_length=$($FASIZE -detailed temp_elongated_subject.fna | cut -f2)
            echo "Modified coords for $subject" >> log.txt
            cat modified_subject_coords >> log.txt

            rm -f temp_elongated_subject.fna
            rm -f temp_standard_subject.fna
            rm -f standard_subject_coords
            rm -f query_fragment.fna
            rm -f query_fragment.faa
            rm -f fragment_align_elong.temp

        else

            echo "No cds found for $subject" >> log.txt
            exit 1

        fi

        # Threshold is the first nucleotide of the codon that gave the amino acid
        # at position sstart in the proteic alignment. Its position is given after
        # elongation of the subject gene

        threshold=$((elongation - 2))

        echo "Variables that are to be sent into query_subject_thresholds : " >> log.txt
        echo "Subject : $subject" >> log.txt
        echo "Threshold : $threshold" >> log.txt
        echo "Elongation : $elongation" >> log.txt
        echo "Sstart_nuc : $sstart_mapped_nuc" >> log.txt
        echo "Qstart : $qstart" >> log.txt
        echo "Elongated length : $elongated_length" >> log.txt

        # The file should remain headerless
        echo -e "$subject\t$threshold\t$elongation\t$sstart_mapped_nuc\t$qstart\t$elongated_length\t$sstart" >> "${query}_subjects_thresholds.tsv"

    done < "query_subjects.tsv"

    
    $FATRANS -stop query.fna query.faa
    $FATRANS elongated_subjects.fna elongated_subjects.faa

    echo -e "qseqid\tqlen_nuc_elong\tsseqid\tslen_nuc_elong\tpident\talign_length\tbar\tgaps\tqstart_nuc_elong\tqend_nuc_elong\tsstart_nuc_elong\tsend_nuc_elong\tevalue\tbitscore\tCIGAR" > "${query}_complete_small_ssearch_elong.tsv"
    echo -e "qseqid\tqlen_nuc_short\tsseqid\tslen_nuc_short\tpident\talign_length\tbar\tgaps\tqstart_nuc_short\tqend_nuc_short\tsstart_nuc_short\tsend_nuc_short\tevalue\tbitscore\tCIGAR" > "${query}_complete_small_ssearch_short.tsv"
    echo -e "qseqid\tqlen_tfastx\tsseqid\tslen_tfastx\tpident\talign_length\tbar\tgaps\tqstart_tfastx\tqend_tfastx\tsstart_tfastx\tsend_tfastx\tevalue\tbitscore\tCIGAR" > "${query}_complete_small_tfastx.tsv"

    # Align the nuc query against the nuc elongated versions of the subjects
    $SSEARCH -3 -n -T8 -XM8G -E 10000 -m8BCL -m "FA ${query}_complete_small_ssearch_elong.aln" query.fna elongated_subjects.fna >> "${query}_complete_small_ssearch_elong.tsv"
    # Align the nuc query against the nuc standard versions of the subjects
    $SSEARCH -3 -n -T8 -XM8G -E 10000 -m8BCL query.fna standard_subjects.fna >> "${query}_complete_small_ssearch_short.tsv"
    # Align the proteic query against the prot elongated versions of the subjects
    $SSEARCH -3 -p -T8 -XM8G -E 10000 -m8BCL -s BL50 -f -11 -g -1 query.faa elongated_subjects.faa >> "${query}_complete_small_tfastx.tsv"
    
    
else 

    touch "${query}_complete_small_ssearch_elong.tsv"
    touch "${query}_complete_small_ssearch_short.tsv"
    touch "${query}_complete_small_tfastx.tsv"
    touch "${query}_subjects_thresholds.tsv"
    touch "${query}_complete_small_ssearch_elong.aln"
    touch "${query}_complete_small_tfastx.aln"

    touch "${query}_fragment_small_ssearch_elong.tsv"
    touch "${query}_fragment_small_ssearch_short.tsv"
    touch "${query}_fragment_small_tfastx.tsv"
    touch "${query}_fragment_small_ssearch_elong.aln"

    echo "No query_subjects.tsv found" >> log.txt


fi