process standardGff {

    containerOptions "--bind ${projectDir}:/DMEL"
    clusterOptions '--partition=common --cpus-per-task=8 --mem=32gb'

    input:
    path gff

    output:
    path "*_standard.gff"


    script:
    """
    #!/usr/bin/bash 

    agat_convert_sp_gxf2gxf.pl --gff /DMEL/input/focal/$gff --out _standard.gff 
    """

}