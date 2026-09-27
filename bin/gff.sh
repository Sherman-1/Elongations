#!/usr/bin/bash

#!/bin/bash

# Function to retrieve unique mRNA IDs and echo them
get_unique_mRNA_ids() {
    # Extract all mRNA IDs from the GFF file and echo them
    awk '$3 == "mRNA" {print $9}' "$1" | grep -ioP 'ID=\K[^;]+' | sort -u
}

# Call the function and capture the result in a variable
mRNA_ids=$(get_unique_mRNA_ids "../input/focal/Dmel.gff")

# Print the unique mRNA IDs
echo "Unique mRNA IDs:"
for id in $mRNA_ids; do
    echo "$id"
done


