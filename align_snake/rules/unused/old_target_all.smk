rule all:
    input:
        expand("{CACHEDIR}/{sample}.h5ad", CACHEDIR=config["CACHEDIR"], sample=SAMPLES),
        # expand("{DATADIR}/align_out/{sample}/Aligned.sortedByCoord.dedup.out.bam.bai", DATADIR=config["DATADIR"], sample=SAMPLES), # umi_tools deduplicated .bam **Note** this is super slow!! Only uncomment if NEEDED
        expand("{DATADIR}/align_out/{sample}/Unmapped_fastqc_out", DATADIR=config["DATADIR"], sample=SAMPLES), #fastQC results for unmapped reads
        expand("{DATADIR}/align_out/{sample}/STARsolo/Unmapped.out.mate2.fastq.gz", DATADIR=config["DATADIR"], sample=SAMPLES), # compress unmapped reads
        expand("{DATADIR}/align_out/{sample}/STARsolo/Solo.out/Gene/filtered/matrix.mtx.gz", DATADIR=config["DATADIR"], sample=SAMPLES), # compressed count matrices
        expand("{DATADIR}/align_out/{sample}/STARsolo/Aligned.sortedByCoord.out.bam.bai", DATADIR=config["DATADIR"], sample=SAMPLES), # index non-dedup .bam

        # expand("{DATADIR}/align_out/{sample}/qualimap_out/qualimapReport.html", DATADIR=config["DATADIR"], sample=SAMPLES), # alignment QC with qualimap)

        expand("{DATADIR}/align_out/{sample}/preTrim_fastqc_R2_out", DATADIR=config["DATADIR"], sample=SAMPLES), # raw R2 fastQC results
        expand("{DATADIR}/align_out/{sample}/postTrim_fastqc_R2_out", DATADIR=config["DATADIR"], sample=SAMPLES), # adapter/polyA/ployG-trimmed R2 fastQC results

        expand("{REFDIR}/{SPECIES}/STAR/SA", REFDIR=REFDIR, SPECIES=SPECIES), # Reference genomes
        expand("{METADIR}/{SRR}.json", METADIR=config["METADIR"], SRR=SRR_LIST), # metadata files
        expand("{METADIR}/merged_metadata.csv", METADIR=config["METADIR"])

        # expand("{DATADIR}/fastqs/{SRR}_{READ}.fastq.gz", DATADIR=config["DATADIR"], SRR=SRR_LIST, READ=["R1", "R2"]) # fastqs
