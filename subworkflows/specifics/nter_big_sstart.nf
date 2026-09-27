include { ssearch_big_sstart as ssearch } from "../../modules/nter/treatBigSstart"
include { parse_big_table } from "../../modules/nter/treatBigSstart"


workflow nter_big_sstart {


    take : 

        big_sstart 
        full_fna
        full_gff

    emit :

        final_nter_big

    main : 

        // TREAT BIG SSTART
        big_tuple = ssearch(big_sstart.flatten(), full_fna, full_gff)
        big_sstart_res = parse_big_table(big_tuple)

        final_nter_big = big_sstart_res.collectFile(name : 'bigFinal', storeDir : 'output/nter/', keepHeader : true, skip : 1)
        
}