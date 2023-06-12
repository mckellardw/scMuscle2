#############################################
## UMI-aware .bam deduplication
#############################################
# Remove reads that don't have a corrected spot/cell barcode with samtools, then remove duplicates w/ **umi-tools**
## Deduplication is done chromosome-by-chromosome to (dramatically) decrease run time/mem usage
rule umitools_dedupBAM:
    input:
        SORTEDBAM = "{DATADIR}/align_out/{sample}/STARsolo/Aligned.sortedByCoord.out.bam"
    output:
        DEDUPBAM = "{DATADIR}/align_out/{sample}/STARsolo/Aligned.sortedByCoord.dedup.out.bam"
    threads:
        config["CORES_LO"]
        #1
    log:
        "{DATADIR}/align_out/{sample}/umitools_dedup.log"
    run:        
        tmp_chemistry = CHEM_DICT[wildcards.sample]
        whitelist = CHEMISTRY_SHEET["whitelist"][tmp_chemistry]

        shell(
            f"""
            bash scripts/split_dedup.sh {input.SORTEDBAM} {whitelist} {threads} {output.DEDUPBAM} {DATADIR}/align_out/{wildcards.sample}/tmp/dedup | tee {log}
            """
        )

# Index the deduplicated .bam file(s)
rule umitools_indexDedupBAM:
    input:
        SORTEDBAM = "{DATADIR}/align_out/{sample}/STARsolo/Aligned.sortedByCoord.dedup.out.bam"
    output:
        BAI = "{DATADIR}/align_out/{sample}/STARsolo/Aligned.sortedByCoord.dedup.out.bam.bai"
    params:
        SAMTOOLS_EXEC = SAMTOOLS_EXEC
    threads:
        config["CORES_LO"]
    shell:
        """
        {params.SAMTOOLS_EXEC} index -@ {threads} {input.SORTEDBAM}
        """
