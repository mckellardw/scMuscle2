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
