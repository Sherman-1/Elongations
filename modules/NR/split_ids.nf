process split_ids {

    input:
    path dataframe
    val n

    output:
    path "*.txt", emit : ids_files
    
    shell:
    '''
    split() {
        if [ $# -ne 2 ]; then
            >&2 echo "Usage: split <arr> <n>"
            return 1
        fi
        local -n arr=$1
        local n=$2
        local len=${#arr[@]}
        local k=$((len / n))
        local m=$((len % n))

        for ((i=0; i<n; i++)); do
            local start=$((i * k + (i<m ? i : m)))
            local end=$(((i+1) * k + (i+1<m ? i+1 : m)))

            for id in "${arr[@]:$start:$((end-start))}"; do
                echo "$id" >> "$i.txt"
		
            done

        done
    }

    readarray -t IDS <<< "$(cut -f 1 !{dataframe} | tail -n +2 | sort -u)"

    split IDS !{n}
    '''
}
