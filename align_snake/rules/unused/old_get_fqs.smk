#############################################
## prefetch (SRA-toolkit)
#############################################
# Download reads as an SRA file into a temporary directory, so they can be formatted as .fastq's
#TODO- add back bam/fastq if statement... :(
rule prefetch:
    output:
        SRA = temp("{TMPDIR}/{SRR}.sra")
    run:
        FORMAT = FORMAT_DICT[wildcards.sample]

        if FORMAT == "fastq":
            shell(
                f"""
                {PREFETCH_EXEC} \
                --verify yes \
                --type sra \
                --max-size 999G \
                --output-file {output.SRA} \
                {wildcards.SRR}
                """
            )
        elif FORMAT == "bam":
            shell(
                f"""
                 {PREFETCH_EXEC} \
                --verify yes \
                --max-size 999999999999 \
                --force ALL \
                --output-directory {TMPDIR}/{wildcards.SRR} \
                --type all \
                --output-file {output.SRA} \
                {wildcards.SRR}
                """
            )
# [SRR for SRR in SRR_LIST if FORMAT_DICT[SRR] in ["fastq","bam"]]

 # Convert SRA files to fastq with fasterq-dump
rule convert_fastqs:
    input:
        SRA = TMPDIR+"/{SRR}.sra"
    output:
        R1_FQ = temp("{DATADIR}/fastqs/{SRR}_1.fastq.gz"),
        R2_FQ = temp("{DATADIR}/fastqs/{SRR}_2.fastq.gz")
    params:
        # TMPDIR = TMPDIR,
        MEMLIMIT = config["MEMLIMIT_LO"],
        # SRA_tmp = "{TMPDIR}/{SRR}.sra",
        # SRA_BAM = "{TMPDIR}/{SRR}.bam",
        OUTDIR = "{DATADIR}/fastqs"
    threads:
        config["CORES_LO"]
    priority:
        42
    run:
        FORMAT = FORMAT_DICT[wildcards.sample]

        if FORMAT == "fastq":
            shell(
                f"""
                {FQD_EXEC} \
                --threads {threads} \
                --mem {params.MEMLIMIT} \
                --outdir {params.OUTDIR} \
                --temp {TMPDIR} \
                --split-files \
                --include-technical \
                {input.SRA}

                pigz -p{threads} --force {params.OUTDIR}/{wildcards.SRR}_*.fastq
                """
            )
        elif FORMAT == "bam":
            shell(
                f"""

                """
            )

# {FFQ} --ncbi {SRR} \
# | xargs curl -o {DATADIR}/fastqs/{SRR}.bam

# Other option:
# ffq --ftp SRR7276476 | grep -Eo '"url": "[^"]*"' | grep -o '"[^"]*"$' | xargs curl -O


# Rule to download fastq files for fastqs hosted on AWS
# rule get_aws_fastqs:
#     output:
#         expand(
#             "{DATADIR}/fastqs/{AWS_FQs}",
#             DATADIR=DATADIR,
#             AWS_FQs = TODO
#         )
#     params:
#         OUTDIR = "{DATADIR}/fastqs/"
#     priority:
#         1
#     run:
#         #TODO

# Rule to download fastq files for fastqs already on the local machine
# rule get_local_fastqs:
#     output:
#         expand(
#             "{DATADIR}/fastqs/{local_FQs}",
#             DATADIR=DATADIR,
#             local_FQs = TODO
#         )
#     params:
#         OUTDIR = "{DATADIR}/fastqs/"
#     priority:
#         1
#     run:
#         TODO

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
