include { small_elongate_and_align } from '../../modules/nter/treatSmallSstart.nf'
include { parse_small_table } from '../../modules/nter/treatSmallSstart.nf'



workflow nter_small_sstart { 


    take : 

        small_sstart 
        full_fna
        full_gff
        index
        
    emit:

        final_nter_small

    main : 


        // TREAT SMALL SSTART
        small_tuple = small_elongate_and_align(small_sstart.flatten(), full_fna, full_gff, index)
        small_sstart_res = parse_small_table(small_tuple)

        final_nter_small = small_sstart_res.collectFile(name : 'smallFinal', storeDir : 'output/nter/', keepHeader : true, skip : 1)

    

        
}