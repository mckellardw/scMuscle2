#############################################
## Trimming and FastQC
#############################################
# Check QC of reads before trimming
rule preTrim_FastQC_R2:
    input:
        MERGED_R2_FQ = "{DATADIR}/align_out/{sample}/tmp/merged_R2.fq.gz"
    output:
        FASTQC_DIR = directory("{DATADIR}/align_out/{sample}/preTrim_fastqc_R2_out"),
        # fastqcReport = ""
    threads:
        config["CORES_LO"]
        # min([config["CORES_LO"],8]) # 8 core max based on recommendations from trim_galore authors
    run:
        shell(
            f"""
            mkdir -p {output.FASTQC_DIR}
            cd {output.FASTQC_DIR}

            {FASTQC_EXEC} \
            --outdir {output.FASTQC_DIR} \
            --threads {threads} \
            {input.MERGED_R2_FQ}
            """
        )

# TSO, polyA, and polyG trimming
rule cutadapt_R2:
    input:
        MERGED_R1_FQ = "{DATADIR}/align_out/{sample}/tmp/merged_R1.fq.gz",
        MERGED_R2_FQ = "{DATADIR}/align_out/{sample}/tmp/merged_R2.fq.gz"
    output:
        FINAL_R1_FQ = "{DATADIR}/align_out/{sample}/tmp/merged_R1_final.fq.gz", #temp()
        FINAL_R2_FQ = "{DATADIR}/align_out/{sample}/tmp/merged_R2_final.fq.gz" #temp()
    params:
        CUTADAPT_EXEC = CUTADAPT_EXEC,
        THREE_PRIME_R2_POLYA = "A"*100, # 100 A-mer
        THREE_PRIME_R2_POLYG = "G"*100, # 100 G-mer
        FIVE_PRIME_R2_TSO = "CCCATGTACTCTGCGTTGATACCACTGCTT", #10x TSO sequence
        FIVE_PRIME_R2_rcTSO = "AAGCAGTGGTATCAACGCAGAGTACATGGG" # rev-comp of 10x TSO sequence
        # FIVE_PRIME_R2 = "TTCGTCACCATAGTTGCGTCTCATGTACCC" #rev 10x TSO sequence
    threads:
        config["CORES_LO"]
        # min([config["CORES_LO"],8]) # 8 core max based on recommendations from trim_galore authors
    log:
        "{DATADIR}/align_out/{sample}/cutadapt.log"
    shell:
        """
        {params.CUTADAPT_EXEC} \
        --minimum-length 18 \
        -A {params.THREE_PRIME_R2_POLYA} \
        -A {params.THREE_PRIME_R2_POLYG} \
        -G {params.FIVE_PRIME_R2_TSO} \
        -G {params.FIVE_PRIME_R2_rcTSO} \
        --pair-filter=any \
        -o {output.FINAL_R1_FQ} \
        -p {output.FINAL_R2_FQ} \
        --cores {threads} \
        {input.MERGED_R1_FQ} {input.MERGED_R2_FQ} 1> {log}
        """

# QC after read trimming
rule postTrim_FastQC_R2:
    input:
        FINAL_R2_FQ = "{DATADIR}/align_out/{sample}/tmp/{sample}_R2_final.fq.gz"
    output:
        FASTQC_DIR = directory("{DATADIR}/align_out/{sample}/postTrim_fastqc_R2_out")
        # fastqcReport = ""
    threads:
        min([config["CORES_LO"],8]) # 8 core max
    run:
        shell(
            f"""
            mkdir -p {output.FASTQC_DIR}
            cd {output.FASTQC_DIR}

            {FASTQC_EXEC} \
            --outdir {output.FASTQC_DIR} \
            --threads {threads} \
            {input.FINAL_R2_FQ}
            """
        )
