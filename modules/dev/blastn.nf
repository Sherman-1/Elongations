process blastn {

    clusterOptions '--partition=common --cpus-per-task=8 --mem=6gb --time=24:00:00'

    label 'blastn'

    input:
    each candidates_and_homologs 
    path full_gff
    path full_fna
    path index

    output:
    path "*_blast.tsv"

    shell:
    '''
    cut -f1,2 !{index} > neighbors.genome

    # Extracting the first column
    query=$(cat !{candidates_and_homologs} | awk '{print $1}')

    second_column=$(cat !{candidates_and_homologs} | awk '{print $2}')

    IFS=';' read -ra homologs <<< ${second_column}

    for homolog in "${homologs[@]}"; do
        
        # Double escaping for nextflow interpreter
        grep -P "\\b$homolog\\b" !{full_gff} | gffread -C - | grep "CDS" | sort -n -k4 -o hom_gff

        num_lines=$(wc -l < hom_gff)

        if [ "$num_lines" -eq 1 ]; then


            bedtools slop -i hom_gff -g neighbors.genome -b 50 | gffread -x - -g !{full_fna} - >> "${query}_subjects.fna"

        else

            rm -f modified_coords
            
            line_number=1
            while read -r line; do

                if [ $line_number -eq 1 ]; then

                    echo "$line" | bedtools slop -i - -g neighbors.genome -l 50 -r 0 >> modified_coords

                elif [ $line_number -eq "$num_lines" ]; then

                    echo "$line" | bedtools slop -i - -g neighbors.genome -l 0 -r 50 >> modified_coords

                else
                    echo "$line" >> modified_coords
                fi

                line_number=$((line_number+1))

                
            done < hom_gff
            
            gffread -x - -g !{full_fna} modified_coords >> "${query}_subjects.fna"

        fi

    done 

    grep -P "\\b$query\\b" !{full_gff} | gffread -C - | grep "CDS" | gffread -x - -g !{full_fna} - >> query.fna

    cat query.fna > test

    echo -e "query\tsubject\tpident\talignlength\tmismatches\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore" >> "${query}_blast.tsv"

    makeblastdb -in "${query}_subjects.fna" -dbtype nucl -out "${query}_db_nucl"
    blastn -task 'dc-megablast' -num_threads 4 -outfmt 6 -db "${query}_db_nucl" -query query.fna -subject_besthit
    '''

}

