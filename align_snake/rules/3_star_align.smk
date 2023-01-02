# Make output directory, align fastqs, and generate raw/filtered feature/cell-barcode matrices
#   Info for STARsolo command line paramaters: https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md
#  *Note:*
#       - `--soloBarcodeReadLength 0` is included to account for files w/ R1 length > than CB+UMI
#       - `priority` should be set higher to ensure
rule STARsolo_align:
    input:
        FINAL_R1_FQ = "{DATADIR}/align_out/{sample}/tmp/merged_R1_final.fq.gz",
        FINAL_R2_FQ = "{DATADIR}/align_out/{sample}/tmp/merged_R2_final.fq.gz",
        REF_LIST = expand("{REFDIR}/{SPECIES}/STAR/Genome", REFDIR=REFDIR, SPECIES=SPECIES) # Reference genomes
    output:
        SORTEDBAM = "{DATADIR}/align_out/{sample}/STARsolo/Aligned.sortedByCoord.out.bam", #TODO: add temp()
        UNMAPPED1 = "{DATADIR}/align_out/{sample}/STARsolo/Unmapped.out.mate1",
        UNMAPPED2 = "{DATADIR}/align_out/{sample}/STARsolo/Unmapped.out.mate2",
        VELDIR = directory("{DATADIR}/align_out/{sample}/STARsolo/Solo.out/Velocyto"),
        GENEDIR = directory("{DATADIR}/align_out/{sample}/STARsolo/Solo.out/Gene"),
        GENEFULLDIR = directory("{DATADIR}/align_out/{sample}/STARsolo/Solo.out/GeneFull"),
        VELMAT = "{DATADIR}/align_out/{sample}/STARsolo/Solo.out/Velocyto/filtered/spliced.mtx",
        GENEMAT = "{DATADIR}/align_out/{sample}/STARsolo/Solo.out/Gene/filtered/matrix.mtx",
        GENEFULLMAT = "{DATADIR}/align_out/{sample}/STARsolo/Solo.out/GeneFull/filtered/matrix.mtx"
    params:
        MEMLIMIT = config["MEMLIMIT_HI"]
    threads:
        config["CORES_MID"]
    priority:
        42
    run:
        tmp_chemistry = CHEM_DICT[wildcards.sample]
        species = SPECIES_DICT[wildcards.sample]

        STAR_REF = f"{REFDIR}/{species}/STAR"
        UMIlen = CHEMISTRY_SHEET["STAR.UMIlen"][tmp_chemistry]
        SOLOtype = CHEMISTRY_SHEET["STAR.soloType"][tmp_chemistry]
        CB_WHITELIST = CHEMISTRY_SHEET["whitelist"][tmp_chemistry]

        print("Using up to " + str(params.MEMLIMIT) + " of memory...")
        shell(
            f"""
            mkdir -p {DATADIR}/align_out/{wildcards.sample}/STARsolo

            {STAR_EXEC} \
            --runThreadN {threads} \
            --outFileNamePrefix {DATADIR}/align_out/{wildcards.sample}/STARsolo/ \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
            --readFilesCommand zcat \
            --genomeDir {STAR_REF} \
            --limitBAMsortRAM={params.MEMLIMIT} \
            --readFilesIn {input.FINAL_R2_FQ} {input.FINAL_R1_FQ} \
            --clipAdapterType CellRanger4 \
            --outReadsUnmapped Fastx \
            --outFilterMultimapNmax 50 \
            --soloUMIlen {UMIlen} \
            --soloType {SOLOtype} \
            --soloCBwhitelist {CB_WHITELIST} \
            --soloBarcodeReadLength 0 \
            --soloCellFilter EmptyDrops_CR \
            --soloFeatures Gene GeneFull Velocyto \
            --soloMultiMappers EM
            """
        )
        # --genomeLoad LoadAndRemove \

# compress outputs from STAR (count matrices)
rule compress_STAR_outs:
    input:
        VELMAT = "{DATADIR}/align_out/{sample}/STARsolo/Solo.out/Velocyto/filtered/spliced.mtx",
        GENEMAT = "{DATADIR}/align_out/{sample}/STARsolo/Solo.out/Gene/filtered/matrix.mtx",
        GENEFULLMAT = "{DATADIR}/align_out/{sample}/STARsolo/Solo.out/GeneFull/filtered/matrix.mtx"
    output:
        VELMAT = "{DATADIR}/align_out/{sample}/STARsolo/Solo.out/Velocyto/filtered/spliced.mtx.gz",
        GENEMAT = "{DATADIR}/align_out/{sample}/STARsolo/Solo.out/Gene/filtered/matrix.mtx.gz",
        GENEFULLMAT = "{DATADIR}/align_out/{sample}/STARsolo/Solo.out/GeneFull/filtered/matrix.mtx.gz"
    params:
        # SOLODIR = "{DATADIR}/align_out/{sample}/STARsolo/Solo.out/"
        VELDIR = directory("{DATADIR}/align_out/{sample}/STARsolo/Solo.out/Velocyto"),
        GENEDIR = directory("{DATADIR}/align_out/{sample}/STARsolo/Solo.out/Gene"),
        GENEFULLDIR = directory("{DATADIR}/align_out/{sample}/STARsolo/Solo.out/GeneFull")
    priority:
        43
    threads:
        1
    run:
        shell(
            f"""
            gzip -qf {params.VELDIR}/*/*.tsv {params.VELDIR}/*/*.mtx
            gzip -qf {params.GENEDIR}/*/*.tsv {params.GENEDIR}/*/*.mtx
            gzip -qf {params.GENEFULLDIR}/*/*.tsv {params.GENEFULLDIR}/*/*.mtx
            """
        )

# Index the .bam produced by STAR
rule indexSortedBAM:
    input:
        SORTEDBAM = "{DATADIR}/align_out/{sample}/STARsolo/Aligned.sortedByCoord.out.bam"
    output:
        BAI = "{DATADIR}/align_out/{sample}/STARsolo/Aligned.sortedByCoord.out.bam.bai"
    params:
        SAMTOOLS_EXEC = SAMTOOLS_EXEC
    threads:
        config["CORES_LO"]
    shell:
        """
        {params.SAMTOOLS_EXEC} index -@ {threads} {input.SORTEDBAM}
        """
