include { small_elongate_and_align } from '../../modules/cter/treatSmallSend.nf'
include { parse_small_table } from '../../modules/cter/treatSmallSend.nf'



workflow cter_small_send { 


    take : 

        small_send 
        full_fna
        full_gff
        index

    emit : 

        final_cter_small

    main : 

        // TREAT SMALL SSTART
        small_tuple = small_elongate_and_align(small_send.flatten(), full_fna, full_gff, index)
        small_send_res = parse_small_table(small_tuple)

        final_cter_small = small_send_res.collectFile(name : 'smallFinal', storeDir : 'output/cter/', keepHeader : true, skip : 1)
}