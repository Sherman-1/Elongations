include { big_diamond as diamond } from '../modules/NR/diamond'
include { parseDiamond } from '../modules/NR/parseDiamond'
include { split_ids } from '../modules/NR/split_ids'
include { split_df } from '../modules/NR/split_df'

workflow NR { 
    
    take:

        focal_proteome
        strain2species_dict
        eukaryotes_dict

    main:

        nr = params.nr

        df = diamond(focal_proteome, nr, 2759)

        split_ids(df, params.n)
	
        subdf = split_df(df, split_ids.out.ids_files.flatten())

        parseDiamond(subdf.flatten(), strain2species_dict, eukaryotes_dict)

        parseDiamond.out.stats.collectFile(name : 'statistics', storeDir : 'test', keepHeader : true, skip : 1)

	    nter_candidates_IDs = parseDiamond.out.nter.collectFile(name : 'nter_candidates', storeDir : 'test')
	    cter_candidates_IDs = parseDiamond.out.cter.collectFile(name : 'cter_candidates', storeDir : 'test')
	

    emit: 

        nter_candidates_IDs
        cter_candidates_IDs

}
