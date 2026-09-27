include { ssearch_big_send as ssearch } from "../../modules/cter/treatBigSend"
include { parse_big_table } from "../../modules/cter/treatBigSend"


workflow cter_big_send {


    take : 

        big_send 
        full_fna
        full_gff

    emit : 

        final_cter_big

    main : 

        // TREAT BIG SSTART
        big_tuple = ssearch(big_send.flatten(), full_fna, full_gff)
        big_send_res = parse_big_table(big_tuple)

        final_cter_big = big_send_res.collectFile(name : 'bigFinal', storeDir : 'output/cter/', keepHeader : true, skip : 1)
        
}