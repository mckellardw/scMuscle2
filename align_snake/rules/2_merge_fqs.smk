# Merge .fastq files (in case more than one sesquencing run was performed)
#     expand("{DATADIR}/fastqs/{SRR}_{READ}.fastq.gz", DATADIR=config["DATADIR"], SRR=SRR_LIST, READ=["R1", "R2"])
rule merge_fastqs:
    input:
       R1_FQ = lambda wildcards: R1_FQs[wildcards.sample],
       R2_FQ = lambda wildcards: R2_FQs[wildcards.sample]
    output:
        MERGED_R1_FQ = temp("{DATADIR}/align_out/{sample}/tmp/{sample}_R1.fq.gz"),
        MERGED_R2_FQ = temp("{DATADIR}/align_out/{sample}/tmp/{sample}_R2.fq.gz")
    params:
        TMP_DIR = "{DATADIR}/align_out/{sample}/tmp"
    priority:
        2
    threads:
        config["CORES_LO"]
    run:
        if len(input.R1_FQ)==1 & len(input.R2_FQ)==1: # shell for single fastq input
            shell("cp {input.R1_FQ} {output.MERGED_R1_FQ}")
            shell("cp {input.R2_FQ} {output.MERGED_R2_FQ}")
        else: # shell enablinging multi-fast input; concatenate inputs
            print("Concatenating",len(input.R1_FQ), ".fastq's for", wildcards.sample)
            shell("mkdir -p {params.TMP_DIR}")
            shell("zcat {input.R1_FQ} > {params.TMP_DIR}/{wildcards.sample}_R1.fq")
            shell("zcat {input.R2_FQ} > {params.TMP_DIR}/{wildcards.sample}_R2.fq")
            shell("pigz -p{threads} {params.TMP_DIR}/{wildcards.sample}_R1.fq {params.TMP_DIR}/{wildcards.sample}_R2.fq")
