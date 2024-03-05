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
        log = "{DATADIR}/align_out/{sample}/fastqc/preTrim_R2.log"
    run:
        shell(
            f"""
            mkdir -p {output.FASTQC_DIR}

            {EXEC['FASTQC']} \
                --outdir {output.FASTQC_DIR} \
                --threads {threads} \
                -a {params.adapters} \
                {input.MERGED_R2_FQ}  \
            1> {log.log}
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
        QUALITY_MIN=20,
        MIN_R2_LENGTH = 12,
        OVERLAP = 5,
        HOMOPOLYMER_ERROR_RATE = 0.2, # default error rate is 0.1
        # THREE_PRIME_R1_POLYA = "A"*100,
        THREE_PRIME_R2_POLYA = "A"*100,
        # THREE_PRIME_R2_POLYG = "G"*100,
        THREE_PRIME_R2_POLYT = "T"*100,
        THREE_PRIME_R2_NEXTERA = "CTGTCTCTTATA", # Nextera sequence
        THREE_PRIME_R2_rcNEXTERA = "TATAAGAGACAG", # Rev Comp of Nextera sequence
        THREE_PRIME_R2_TSO = "AAGCTGGTATCAACGCAGAGTGAATGGG", # SlideSeq TSO - remove any polyadenylated TSOs
        THREE_PRIME_R2_TXG_TSO = "AAGCAGTGGTATCAACGCAGAGTACATGGG", # 10x TSO - remove any polyadenylated TSOs
        THREE_PRIME_R2_ILLUMINA_UNI = "AGATCGGAAGAG", # Illumina Universal
        FIVE_PRIME_R2_TSO = "CCCATTCACTCTGCGTTGATACCAGCTT", # rev comp of SlideSeq TSO
        FIVE_PRIME_R2_TXG_TSO = "CCCATGTACTCTGCGTTGATACCACTGCTT", # rev-comp of 10x TSO sequence
    threads:
        config["CORES_LO"]
    log:
        log = "{DATADIR}/align_out/{sample}/cutadapt.log"
    run:
        shell(
            f"""
            {EXEC['CUTADAPT']} \
                --minimum-length {params.MIN_R2_LENGTH} \
                --quality-cutoff {params.QUALITY_MIN} \
                --overlap {params.OVERLAP} \
                --match-read-wildcards \
                --nextseq-trim=20 \
                -A "{params.THREE_PRIME_R2_POLYA};max_error_rate={params.HOMOPOLYMER_ERROR_RATE}" \
                -A "{params.THREE_PRIME_R2_POLYT};max_error_rate={params.HOMOPOLYMER_ERROR_RATE}" \
                -A {params.THREE_PRIME_R2_TSO} \
                -A {params.THREE_PRIME_R2_TXG_TSO} \
                -A {params.THREE_PRIME_R2_NEXTERA} \
                -A {params.THREE_PRIME_R2_rcNEXTERA} \
                -A {params.THREE_PRIME_R2_ILLUMINA_UNI} \
                -G {params.FIVE_PRIME_R2_TSO} \
                -G {params.FIVE_PRIME_R2_TXG_TSO} \
                --pair-filter=any \
                -o {output.FINAL_R1_FQ} \
                -p {output.FINAL_R2_FQ} \
                --cores {threads} \
                {input.MERGED_R1_FQ} {input.MERGED_R2_FQ} \
            1> {log.log}
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
        log = "{DATADIR}/align_out/{sample}/fastqc/postTrim_R2.log"
    run:
        shell(
            f"""
            mkdir -p {output.FASTQC_DIR}

            {EXEC['FASTQC']} \
                --outdir {output.FASTQC_DIR} \
                --threads {threads} \
                -a {params.adapters} \
                {input.FINAL_R2_FQ} \
            1> {log.log}
            """
        )
