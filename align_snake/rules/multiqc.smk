#############################################
## Aggregate QC information with multiQC
#############################################

rule multiqc_final:
    input:
        expand()
    output:
        MQC_DIR = directory("{DATADIR}/multiqc_out/")
    threads:
        config["CORES_HI"]
    run:
