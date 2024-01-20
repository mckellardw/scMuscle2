#############################################
## prefetch (SRA-toolkit)
#############################################
# Download reads as an SRA file into a temporary directory, so they can be formatted as .fastq's

 # Download and convert .sra files to .fastq with fasterq-dump
 #TODO- replace SRA nonsense with [fastq-dl](https://github.com/rpetit3/fastq-dl)
rule get_fastqs:
    input:
        SRR_LIST = "{DATADIR}/align_out/{sample}/SRR_list.txt"
    output:
        MERGED_R1_FQ = "{DATADIR}/align_out/{sample}/tmp/merged_R1.fq.gz", #temp()
        MERGED_R2_FQ = "{DATADIR}/align_out/{sample}/tmp/merged_R2.fq.gz"
    params:
        MEMLIMIT = config["MEMLIMIT_LO"]
    threads:
        config["CORES_LO"]
    priority:
        42
    run:
        FORMAT    = FORMAT_DICT[wildcards.sample]
        CHEMISTRY = CHEM_DICT[wildcards.sample]
        # [SRR for SRR in SRR_LIST if FORMAT_DICT[SRR] in ["fastq","bam"]]

        if FORMAT == "fastq":
            shell(
                f"""
                cd {DATADIR}/align_out/{wildcards.sample}
                python {PRODIR}/align_snake/scripts/get_sra_fq.py {EXEC['PREFETCH']} {EXEC['FASTERQDUMP']} {threads} {params.MEMLIMIT}
                """
            )
        elif FORMAT == "bam":
            # pull .bam file
            shell(
                f"""
                cd {DATADIR}/align_out/{wildcards.sample}
                bash {PRODIR}/align_snake/scripts/get_sra_bam.sh {EXEC['PREFETCH']} {EXEC['BAM2FQ']} {threads} {wildcards.sample}
                """
            )
        elif FORMAT == "aws": 
            LINKS = LINK_DICT[wildcards.sample]
            LINKS = LINKS.replace(";", " " ) #make sure file list is space-delimited
            shell(
                f"""
                cd {DATADIR}/align_out/{wildcards.sample}
                mkdir -p tmp/pulled_fqs
                cd tmp/pulled_fqs
                rm ./*
                wget {LINKS}
                zcat ./*R1*.f*.gz > ../merged_R1.fq
                zcat ./*R2*.f*.gz > ../merged_R2.fq
                pigz -p{threads} ../merged_R*.fq
                rm ./*
                """
            )
        elif FORMAT == "local": #TODO
            R1_LIST = str()
            R2_LIST = str()
            shell(
                f"""
                cd {DATADIR}/align_out/{wildcards.sample}
                zcat {R1_LIST} > tmp/merged_R1.fq
                zcat {R2_LIST} > tmp/merged_R2.fq
                pigz -p{threads} tmp/*.fq
                """
            )
        elif FORMAT == "ngdc_fastq": #TODO
            LINKS = LINK_DICT[wildcards.sample]
            LINKS = LINKS.replace(";", " " ) #make sure file list is space-delimited
            shell(
                f"""
                cd {DATADIR}/align_out/{wildcards.sample}
                mkdir -p tmp/pulled_fqs
                cd tmp/pulled_fqs
                rm ./*
                wget {LINKS}
                zcat ./*R1*.f*.gz > ../merged_R1.fq
                zcat ./*R2*.f*.gz > ../merged_R2.fq
                pigz -p{threads} ../merged_R*.fq
                rm ./*
                """
            )

# {FFQ} --ncbi {SRR} \
# | xargs curl -o {DATADIR}/fastqs/{SRR}.bam

# Other option:
# ffq --ftp SRR7276476 | grep -Eo '"url": "[^"]*"' | grep -o '"[^"]*"$' | xargs curl -O
