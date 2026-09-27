// Includes subworkflows
include { nter_big_sstart } from './specifics/nter_big_sstart'
include { nter_small_sstart } from './specifics/nter_small_sstart' 

include { cter_big_send } from './specifics/cter_big_send'
include { cter_small_send } from './specifics/cter_small_send'

// Include modules
include { seqkitGrep as grep_nter } from '../modules/grepSeq'
include { seqkitGrep as grep_cter } from '../modules/grepSeq'
include { nter as parseNter } from "../modules/parseNeighborsBlastp"
include { cter as parseCter } from "../modules/parseNeighborsBlastp"
include { blastp_strict as blastp_strict_nter } from "../modules/blastp"
include { blastp_strict as blastp_strict_cter } from "../modules/blastp"
include { ssearch as ssearch_nter } from "../modules/blastp"
include { ssearch as ssearch_cter } from "../modules/blastp"



workflow search_homologs { 

    take: 

        nter_candidates_IDs
        cter_candidates_IDs 
        local_db
        focal_proteome
        full_gff
        full_fna
        index
        geneCoords
        local_proteome
        

    main:


        // #######################  
        // #                     #
        // #   Nter candidates   #
        // #                     #
        // #######################

        // Cherche les séquences
        // Pardon je dois relancer ... 
        nter_candidates_proteins = grep_nter(nter_candidates_IDs, focal_proteome)
        nter_candidates_homologs = blastp_strict_nter(nter_candidates_proteins, local_db)
        
        ( small_sstart, big_sstart, bad_species, full_blast ) = parseNter(nter_candidates_homologs, geneCoords)

        // Send queries and their homologs to the corresponding subworkflow
        final_nter_small = nter_small_sstart(small_sstart, full_fna, full_gff, index)
        final_nter_big = nter_big_sstart(big_sstart, full_fna, full_gff)


        // Artefact, keep as is for now
        //nter_candidates_homologs = ssearch_nter(nter_candidates_proteins, local_proteome, 75, 200)


        // #######################
        // #                     #
        // #   Cter candidates   #
        // #                     #
        // #######################
        
        // Cherche les séquences
        cter_candidates_proteins = grep_cter(cter_candidates_IDs, focal_proteome)
        cter_candidates_homologs = blastp_strict_cter(cter_candidates_proteins, local_db)
        
        ( small_send, big_send, bad_species, full_blast ) = parseCter(cter_candidates_homologs, geneCoords)

        final_cter_small = cter_small_send(small_send, full_fna, full_gff, index)
        final_cter_big = cter_big_send(big_send,full_fna, full_gff )


        // Artefact, keep as is for now
        //cter_candidates_homologs = ssearch_cter(cter_candidates_proteins, local_proteome, 80, 100)
        




        

    
        

}
