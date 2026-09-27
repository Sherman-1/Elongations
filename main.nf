include { NR } from "./subworkflows/A2_NR" 
include { search_homologs } from "./subworkflows/A3_search_extensions"
include { create_taxon_maps } from "./modules/create_taxon_maps"
include { preprocess_input_data } from "./subworkflows/A1_preprocess"


workflow {
   // Preprocess entry data, extract proteomes and create local database
   (full_gff, full_fna, index, focal_proteome, local_db, geneCoords, local_proteome) = preprocess_input_data()

   (strain2species, eukaryotes) = create_taxon_maps()

   (nter_candidates_IDs, cter_candidates_IDs) = NR(focal_proteome, strain2species, eukaryotes)

   search_homologs(nter_candidates_IDs, cter_candidates_IDs, local_db, focal_proteome, full_gff, full_fna, index, geneCoords, local_proteome)

}
