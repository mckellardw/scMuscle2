#############################################
## Ambient RNA decontamination
#############################################
rule ambient_rna_decon:
    input:
        GENEFULL_FILT_MAT = "{DATADIR}/align_out/{sample}/STARsolo/Solo.out/GeneFull/filtered/matrix.mtx.gz",
        FILT_CELLS = "{DATADIR}/align_out/{sample}/STARsolo/Solo.out/GeneFull/filtered/barcodes.tsv.gz",
        FILT_FEATS = "{DATADIR}/align_out/{sample}/STARsolo/Solo.out/GeneFull/filtered/features.tsv.gz",
        GENEFULL_RAW_MAT = "{DATADIR}/align_out/{sample}/STARsolo/Solo.out/GeneFull/raw/matrix.mtx.gz",
        RAW_CELLS = "{DATADIR}/align_out/{sample}/STARsolo/Solo.out/GeneFull/raw/barcodes.tsv.gz",
        RAW_FEATS = "{DATADIR}/align_out/{sample}/STARsolo/Solo.out/GeneFull/raw/features.tsv.gz"
    output:
        SOUPX_MAT = "{DATADIR}/align_out/{sample}/STARsolo/Solo.out/GeneFull/soupx/matrix.mtx.gz",
        SOUPX_CELLS = "{DATADIR}/align_out/{sample}/STARsolo/Solo.out/GeneFull/raw/barcodes.tsv.gz",
        SOUPX_FEATS = "{DATADIR}/align_out/{sample}/STARsolo/Solo.out/GeneFull/raw/features.tsv.gz"
    params:
    threads:
        config["CORES_MID"]
    script:
        """
        Rscript scripts/starsolo_soupx.R {DATADIR}/align_out/{wildcards.sample}/STARsolo/Solo.out/GeneFull/ {threads}
        """
