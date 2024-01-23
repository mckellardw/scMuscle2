#############################################
## Trimming and FastQC
#############################################
# Check QC of reads before trimming
rule preTrim_FastQC_R2:
    input:
        MERGED_R2_FQ = "{DATADIR}/align_out/{sample}/tmp/merged_R2.fq.gz"
    output:
        FASTQC_DIR = directory("{DATADIR}/align_out/{sample}/fastqc/preTrim_R2"),
    params:
        adapters = config['FASTQC_ADAPTERS']
    threads:
        config["CORES_LO"]
    log:
        "{DATADIR}/align_out/{sample}/fastqc/preTrim_R2.log"
    run:
        shell(
            f"""
            mkdir -p {output.FASTQC_DIR}

            {EXEC['FASTQC']} \
            --outdir {output.FASTQC_DIR} \
            --threads {threads} \
            -a {params.adapters} \
            {input.MERGED_R2_FQ}
            """
        )

# TSO, polyA, and polyG trimming
rule cutadapt:
    input:
        MERGED_R1_FQ = "{DATADIR}/align_out/{sample}/tmp/merged_R1.fq.gz",
        MERGED_R2_FQ = "{DATADIR}/align_out/{sample}/tmp/merged_R2.fq.gz"
    output:
        FINAL_R1_FQ = "{DATADIR}/align_out/{sample}/tmp/merged_R1_final.fq.gz", #temp()
        FINAL_R2_FQ = "{DATADIR}/align_out/{sample}/tmp/merged_R2_final.fq.gz" #temp()
    params:
        MIN_R2_LENGTH = 16,
        THREE_PRIME_R2_POLYA = "A"*100, # 100 A-mer
        THREE_PRIME_R2_POLYG = "G"*100, # 100 G-mer
        FIVE_PRIME_R2_TSO = "CCCATGTACTCTGCGTTGATACCACTGCTT", #10x TSO sequence
        FIVE_PRIME_R2_rcTSO = "AAGCAGTGGTATCAACGCAGAGTACATGGG" # rev-comp of 10x TSO sequence
        # FIVE_PRIME_R2 = "TTCGTCACCATAGTTGCGTCTCATGTACCC" #rev 10x TSO sequence
    threads:
        config["CORES_LO"]
    log:
        "{DATADIR}/align_out/{sample}/cutadapt.log"
    run:
        shell(
            f"""
            {EXEC['CUTADAPT']} \
            --minimum-length {PARAMS.MIN_R2_LENGTH} \
            -A {params.THREE_PRIME_R2_POLYA} \
            -A {params.THREE_PRIME_R2_POLYG} \
            -G {params.FIVE_PRIME_R2_TSO} \
            -G {params.FIVE_PRIME_R2_rcTSO} \
            --pair-filter=any \
            -o {output.FINAL_R1_FQ} \
            -p {output.FINAL_R2_FQ} \
            --cores {threads} \
            {input.MERGED_R1_FQ} {input.MERGED_R2_FQ} \
            1> {log}
            """
        )

# QC after read trimming
rule postTrim_FastQC_R2:
    input:
        FINAL_R2_FQ = "{DATADIR}/align_out/{sample}/tmp/merged_R2_final.fq.gz"
    output:
        FASTQC_DIR = directory("{DATADIR}/align_out/{sample}/fastqc/postTrim_R2")
    params:
        adapters = config['FASTQC_ADAPTERS']
    threads:
        min([config["CORES_LO"], 8]) # 8 core max
    log:
        "{DATADIR}/align_out/{sample}/fastqc/postTrim_R2.log"
    run:
        shell(
            f"""
            mkdir -p {output.FASTQC_DIR}

            {EXEC['FASTQC']} \
            --outdir {output.FASTQC_DIR} \
            --threads {threads} \
            -a {params.adapters} \
            {input.FINAL_R2_FQ}
            """
        )
