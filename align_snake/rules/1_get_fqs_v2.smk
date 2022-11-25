#############################################
## prefetch (SRA-toolkit)
#############################################
# Download reads as an SRA file into a temporary directory, so they can be formatted as .fastq's

 # Convert SRA files to fastq with fasterq-dump
rule get_fastqs:
    input:
        SRR_LIST = "{DATADIR}/align_out/{sample}/SRR_list.txt"
    output:
        MERGED_R1_FQ = temp("{DATADIR}/align_out/{sample}/tmp/merged_R1.fq.gz"),
        MERGED_R2_FQ = temp("{DATADIR}/align_out/{sample}/tmp/merged_R2.fq.gz")
    params:
        MEMLIMIT = config["MEMLIMIT_LO"]
    threads:
        config["CORES_LO"]
    priority:
        42
    run:
        FORMAT = FORMAT_DICT[wildcards.sample]
        CHEMISTRY = CHEM_DICT[wildcards.sample]
        # [SRR for SRR in SRR_LIST if FORMAT_DICT[SRR] in ["fastq","bam"]]

        if FORMAT == "fastq":
            shell(
            f"""
            cd {DATADIR}/align_out/{wildcards.sample}
            python {PRODIR}/align_snake/scripts/get_sra_fq.py {PREFETCH_EXEC} {FQD_EXEC} {threads} {params.MEMLIMIT}
            """
            )
        elif FORMAT == "bam":
            # pull .bam file
            shell(
                f"""
                cd {DATADIR}/align_out/{wildcards.sample}
                bash {PRODIR}/align_snake/scripts/get_sra_bam.sh {PREFETCH_EXEC} {BAM2FQ_EXEC} {threads} {wildcards.sample}
                """
            )
        elif FORMAT == "aws":
            print("TODO")
            # download with wget
            # merge into GSM####.fastq
            # compress
        elif FORMAT == "local":
            print("TODO")
            # merge into GSM####.fastq
            # compress

# {FFQ} --ncbi {SRR} \
# | xargs curl -o {DATADIR}/fastqs/{SRR}.bam

# Other option:
# ffq --ftp SRR7276476 | grep -Eo '"url": "[^"]*"' | grep -o '"[^"]*"$' | xargs curl -O
