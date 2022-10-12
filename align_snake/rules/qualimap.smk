#############################################
## QC on STAR outputs
#############################################

## qualimap on aligned reads
rule qualimapQC:
    input:
        SORTEDBAM = "{DATADIR}/align_out/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        qualimapDir = directory("{DATADIR}/align_out/{sample}/qualimap_out"),
        fastqcReport = "{DATADIR}/align_out/{sample}/qualimap_out/qualimapReport.html"
    params:
        QUALIMAP_EXEC = QUALIMAP_EXEC,
        MEMLIMIT = config["MEMLIMIT_LO"]
        # GENES_GTF = config["GENES_GTF"] #TODO Pull from downloaded gget files
    threads:
        1
        # config["CORES_LO"]
    shell:
        """
        mkdir -p {output.qualimapDir}
        cd {output.qualimapDir}

        {params.QUALIMAP_EXEC} rnaseq \
        -bam {input.SORTEDBAM} \
        -gtf {params.GENES_GTF} \
        --sequencing-protocol strand-specific-forward \
        --sorted \
        --java-mem-size={params.MEMLIMIT} \
        -DATADIR {output.qualimapDir} \
        -outformat html
        """
