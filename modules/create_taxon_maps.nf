process create_taxon_maps {

    clusterOptions '--partition=common --cpus-per-task=8'


    label 'create_taxon_maps'

    output:
    path "strain2species.csv"
    path "eukaryotes.csv"

    script:
    """
    TAXONKIT=${params.taxonkit}

    # List of eukaryotes taxids 
    \$TAXONKIT list --ids 2759 --indent "" | awk '\$1 != "" {print \$1";1"}' > eukaryotes.csv

    \$TAXONKIT list --ids 2759 --indent "" | \$TAXONKIT filter -L species --discard-noranks \\
        | \$TAXONKIT reformat -I 1 --format "{s}" --show-lineage-taxids \\
        | cut -f 1,3 | awk '{print \$1";"\$2}' > strain2species.csv

    # Species to themselves in the same file
    \$TAXONKIT list --ids 2759 --indent "" | \$TAXONKIT filter -E species \\
        | awk '{print \$1";"\$1}' >> strain2species.csv

    """

}