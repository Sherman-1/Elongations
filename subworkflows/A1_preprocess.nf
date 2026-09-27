include { extract_protein as extract_focal_proteome } from '../modules/preprocess/extract_protein'
include { extract_protein as extract_neighbors_proteome } from '../modules/preprocess/extract_protein'
include { diamond_makedb } from '../modules/preprocess/makedb'
include { blast_makedb } from '../modules/preprocess/makedb'
include { faidx } from '../modules/preprocess/faidx'
include { concat as concat_fasta } from '../modules/preprocess/concat'
include { concat as concat_gff } from '../modules/preprocess/concat'
include { concat as concat_prot } from '../modules/preprocess/concat'
include { sample_proteins } from '../modules/preprocess/sample_proteins'
include { getGeneCoords } from '../modules/preprocess/getGeneCoords'
include { checkAnnot }  from '../modules/preprocess/check_annotation'


workflow preprocess_input_data { 

    main:
    
        // Collect genomic annotations and sequences for all species
        neighbors_paired_annot = Channel.fromFilePairs("input/neighbors/*{.gff,.fna}", flat : true)
        focal_paired_annot = Channel.fromFilePairs("input/focal/*{.gff,.fna}", flat : true)

        // Check if the input files are correctly formatted
        //checkAnnot()

        // Collect all gff and fna files separately for ulterior use
        full_gff = concat_gff(Channel.fromPath("input/*/*.gff").collect())
        full_fna = concat_fasta(Channel.fromPath("input/*/*.fna").collect())
        index = faidx(full_fna)
        geneCoords = getGeneCoords(full_gff)

        // Sample or extract proteome of focal species
        focal_proteins = extract_focal_proteome(focal_paired_annot)

        //focal_proteins = sample_proteins(focal_paired_annot, 2000)

        // Keep proteome of all species and concat into one file
        neighbors_flat_proteins = extract_neighbors_proteome(neighbors_paired_annot)
        neighbors_concatenated_proteins = concat_prot(neighbors_flat_proteins.collect())

        // Create a local database from the full proteome
        local_db = blast_makedb(neighbors_concatenated_proteins)

    emit:

        full_gff
        full_fna
        index
        focal_proteins
        local_db
        geneCoords
        neighbors_concatenated_proteins


}
