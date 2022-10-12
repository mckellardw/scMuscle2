#############################################
## UMI-aware .bam deduplication
#############################################
# Remove reads that don't have a corrected spot/cell barcode with samtools, then remove duplicates w/ **umi-tools**
## High mem usage? Check here! https://umi-tools.readthedocs.io/en/latest/faq.html
## **WARNING** this step is suuuupppppeeerrrr sloowwww. Don't run it if you don't need to!
rule umitools_dedupBAM:
    input:
        CB_WHITELIST = CHEMISTRY_SHEET, #TODO: fix this...
        SORTEDBAM = "{DATADIR}/align_out/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        DEDUPBAM = "{DATADIR}/align_out/{sample}/Aligned.sortedByCoord.dedup.out.bam",
        TMPBAM = temp("{DATADIR}/align_out/{sample}/tmp.bam")
    params:
        SAMTOOLS_EXEC = SAMTOOLS_EXEC,
        UMITOOLS_EXEC = UMITOOLS_EXEC,
        OUTPUT_PREFIX="{DATADIR}/align_out/{sample}/umitools_dedup/{sample}"
        # TMPBAM = "{DATADIR}/align_out/{sample}/tmp.bam"
    threads:
        config["CORES_LO"]
        #1
    log:
        "{DATADIR}/align_out/{sample}/umitools_dedup/dedup.log"
    shell:
        """
        TODO
        """

rule umitools_indexDedupBAM:
    input:
        SORTEDBAM = "{DATADIR}/align_out/{sample}/Aligned.sortedByCoord.dedup.out.bam"
    output:
        BAI = "{DATADIR}/align_out/{sample}/Aligned.sortedByCoord.dedup.out.bam.bai"
    params:
        SAMTOOLS_EXEC = SAMTOOLS_EXEC
    threads:
        config["CORES_LO"]
    shell:
        """
        {params.SAMTOOLS_EXEC} index -@ {threads} {input.SORTEDBAM}
        """
