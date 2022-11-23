#############################################
## prefetch (SRA-toolkit)
#############################################
# Download reads as an SRA file into a temporary directory, so they can be formatted as .fastq's

 # Convert SRA files to fastq with fasterq-dump
rule get_fastqs:
    input:
        SRR_LIST = "{DATADIR}/align_out/{sample}/SRR_list.txt"
    output:
        MERGED_R1_FQ = temp("{DATADIR}/align_out/{sample}/tmp/{sample}_R1.fq.gz"),
        MERGED_R2_FQ = temp("{DATADIR}/align_out/{sample}/tmp/{sample}_R2.fq.gz")
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
            """
            cd {DATADIR}/align_out/{sample}
            python scripts/get_sra_fq.py {PREFETCH_EXEC} {FQD_EXEC} {threads} {params.MEMLIMT}
            """
            )
        elif FORMAT == "bam":
            # pull .bam file
            shell(
                f"""
                bash scripts/get_sra_bam.sh {PREFETCH_EXEC} {BAM2FQ_EXEC} {threads} {wildcards.sample}
                """
            )
        elif FORMAT == "aws":
            # download with wget
            # merge into GSM####.fastq
            # compress
        elif FORMAT == "local":
            # merge into GSM####.fastq
            # compress

# {FFQ} --ncbi {SRR} \
# | xargs curl -o {DATADIR}/fastqs/{SRR}.bam

# Other option:
# ffq --ftp SRR7276476 | grep -Eo '"url": "[^"]*"' | grep -o '"[^"]*"$' | xargs curl -O


# Check fastq sizes, then rename/label files for R1 & R2 (discard I1)
#TODO- account for dual index (v3.1) SRRs!!
#TODO- account for paired-end samples (same R1 and R2 lengths)
rule rename_SRR_fastqs:
    input:
        R1_FQ = "{DATADIR}/fastqs/{SRR}_1.fastq.gz",
        R2_FQ = "{DATADIR}/fastqs/{SRR}_2.fastq.gz"
    output:
        R1_FQ = temp("{DATADIR}/fastqs/{SRR}_R1.fastq.gz"),
        R2_FQ = temp("{DATADIR}/fastqs/{SRR}_R2.fastq.gz")
    params:
        OUTDIR = "{DATADIR}/fastqs/"
    priority:
        1
    run:
        #TODO- rewrite this into a function that doesn't depend on file size, and can actually do this programmatically...
        # if params.UPLOADTYPE == "fastq":
            # Check for file sizes (R2 is biggest, then R1, then I2 & I1)
            if path.exists(f"{DATADIR}/fastqs/{wildcards.SRR}_4.fastq.gz"):
                sra_sizes = {
                    f"{DATADIR}/fastqs/{wildcards.SRR}_1.fastq.gz": path.getsize(f"{DATADIR}/fastqs/{wildcards.SRR}_1.fastq.gz"),
                    f"{DATADIR}/fastqs/{wildcards.SRR}_2.fastq.gz": path.getsize(f"{DATADIR}/fastqs/{wildcards.SRR}_2.fastq.gz"),
                    f"{DATADIR}/fastqs/{wildcards.SRR}_3.fastq.gz": path.getsize(f"{DATADIR}/fastqs/{wildcards.SRR}_3.fastq.gz"),
                    f"{DATADIR}/fastqs/{wildcards.SRR}_4.fastq.gz": path.getsize(f"{DATADIR}/fastqs/{wildcards.SRR}_4.fastq.gz")
                }
                rename_dict = {
                    "R2": max(sra_sizes),
                    "R1": list(sra_sizes.keys())[list(sra_sizes.values()).index(sorted(sra_sizes.values())[-3])],
                    "I2": list(sra_sizes.keys())[list(sra_sizes.values()).index(sorted(sra_sizes.values())[-2])],
                    "I1": min(sra_sizes)
                }

                remove(rename_dict["I2"])
                remove(rename_dict["I1"])
                rename(rename_dict["R1"], f"{DATADIR}/fastqs/{wildcards.SRR}_R1.fastq.gz")
                rename(rename_dict["R2"], f"{DATADIR}/fastqs/{wildcards.SRR}_R2.fastq.gz")
            elif path.exists(f"{DATADIR}/fastqs/{wildcards.SRR}_3.fastq.gz"):
                sra_sizes = {
                    f"{DATADIR}/fastqs/{wildcards.SRR}_1.fastq.gz": path.getsize(f"{DATADIR}/fastqs/{wildcards.SRR}_1.fastq.gz"),
                    f"{DATADIR}/fastqs/{wildcards.SRR}_2.fastq.gz": path.getsize(f"{DATADIR}/fastqs/{wildcards.SRR}_2.fastq.gz"),
                    f"{DATADIR}/fastqs/{wildcards.SRR}_3.fastq.gz": path.getsize(f"{DATADIR}/fastqs/{wildcards.SRR}_3.fastq.gz")
                }
                rename_dict = {
                    "R2": max(sra_sizes),
                    "R1": list(sra_sizes.keys())[list(sra_sizes.values()).index(sorted(sra_sizes.values())[-2])],
                    "I1": min(sra_sizes)
                }

                remove(rename_dict["I1"])
                rename(rename_dict["R1"], f"{DATADIR}/fastqs/{wildcards.SRR}_R1.fastq.gz")
                rename(rename_dict["R2"], f"{DATADIR}/fastqs/{wildcards.SRR}_R2.fastq.gz")
            else:
                sra_sizes = {
                    f"{DATADIR}/fastqs/{wildcards.SRR}_1.fastq.gz": path.getsize(f"{DATADIR}/fastqs/{wildcards.SRR}_1.fastq.gz"),
                    f"{DATADIR}/fastqs/{wildcards.SRR}_2.fastq.gz": path.getsize(f"{DATADIR}/fastqs/{wildcards.SRR}_2.fastq.gz")
                }

                rename_dict = {
                    "R2": max(sra_sizes),
                    "R1": min(sra_sizes)
                }

                rename(rename_dict["R1"], f"{DATADIR}/fastqs/{wildcards.SRR}_R1.fastq.gz")
                rename(rename_dict["R2"], f"{DATADIR}/fastqs/{wildcards.SRR}_R2.fastq.gz")
        # elif params.UPLOADTYPE == "bam":
        #     #rename R1 & R2 files (can use R1 and R2 from file names here, not dependent on file sizes)
        #     shell(
        #         """
        #         mv {input.R1_FQ} {output.R1_FQ}
        #         mv {input.R2_FQ} {output.R2_FQ}
        #         """
        #     )
