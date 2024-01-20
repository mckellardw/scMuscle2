#############################################
## Unmapped read analyses
#############################################

# Compress unmapped reads; switch names and add `.fastq` file extension because of STAR weirdness
rule unmapped_compress:
    input:
        UNMAPPED1 = "{DATADIR}/align_out/{sample}/STARsolo/Unmapped.out.mate1",
        UNMAPPED2 = "{DATADIR}/align_out/{sample}/STARsolo/Unmapped.out.mate2"
    output:
        UNMAPPED1_FQ = "{DATADIR}/align_out/{sample}/STARsolo/Unmapped.out.mate1.fastq.gz",
        UNMAPPED2_FQ = "{DATADIR}/align_out/{sample}/STARsolo/Unmapped.out.mate2.fastq.gz"
    threads:
        config["CORES_LO"]
    run:
        shell(
            f"""
            mv {input.UNMAPPED1} {input.UNMAPPED2}.fastq
            mv {input.UNMAPPED2} {input.UNMAPPED1}.fastq

            {EXEC['PIGZ']} -p {threads} {input.UNMAPPED1}.fastq {input.UNMAPPED2}.fastq
            """
        )

# Run fastqc on unmapped reads
rule unmapped_fastqc:
    input:
        UNMAPPED1_FQ = "{DATADIR}/align_out/{sample}/STARsolo/Unmapped.out.mate1.fastq.gz",
        UNMAPPED2_FQ = "{DATADIR}/align_out/{sample}/STARsolo/Unmapped.out.mate2.fastq.gz"
    output:
        FASTQC_DIR = directory("{DATADIR}/align_out/{sample}/fastqc/unmapped_R2")
    threads:
        config["CORES_LO"]
    run:
        shell(
            f"""
            mkdir -p {output.FASTQC_DIR}

            {EXEC['FASTQC']} \
            -o {output.FQC_DIR} \
            -t {threads} \
            {input.UNMAPPED2_FQ}
            """
        )
