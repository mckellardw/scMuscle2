#############################################
## Ambient RNA decontamination
#############################################
rule ambient_rna_decon_soupx:
    input:
        GENEFULL_FILT_MAT = DATADIR+"/align_out/{sample}/STARsolo/Solo.out/GeneFull/filtered/matrix.mtx.gz",
        # FILT_CELLS = "{DATADIR}/align_out/{sample}/STARsolo/Solo.out/GeneFull/filtered/barcodes.tsv.gz",
        # FILT_FEATS = "{DATADIR}/align_out/{sample}/STARsolo/Solo.out/GeneFull/filtered/features.tsv.gz",
        GENEFULL_RAW_MAT = DATADIR+"/align_out/{sample}/STARsolo/Solo.out/GeneFull/raw/matrix.mtx.gz"
        # RAW_CELLS = "{DATADIR}/align_out/{sample}/STARsolo/Solo.out/GeneFull/raw/barcodes.tsv.gz",
        # RAW_FEATS = "{DATADIR}/align_out/{sample}/STARsolo/Solo.out/GeneFull/raw/features.tsv.gz"
    output:
        SOUPX_MAT = DATADIR+"/align_out/{sample}/STARsolo/Solo.out/GeneFull/soupx/matrix.mtx.gz",
        SOUPX_CELLS = DATADIR+"/align_out/{sample}/STARsolo/Solo.out/GeneFull/soupx/barcodes.tsv.gz",
        SOUPX_FEATS = DATADIR+"/align_out/{sample}/STARsolo/Solo.out/GeneFull/soupx/features.tsv.gz"
    params:
    threads:
        config["CORES_MID"]
    run:
        shell(
            f"""
            Rscript scripts/starsolo_soupx.R {DATADIR}/align_out/{wildcards.sample}/STARsolo/Solo.out/GeneFull/ {threads}
            """
        )

rule cache_soupx_h5ad:
    input:
        GENEFULLMAT = DATADIR+"/align_out/{sample}/STARsolo/Solo.out/GeneFull/soupx/matrix.mtx.gz"
    output:
        H5AD = DATADIR+"/align_out/{sample}/STARsolo/Solo.out/GeneFull/soupx/deconQC.h5ad"
    threads:
        1
    run:
        shell(
            f"""
            python scripts/cache_h5ad.py {DATADIR}/align_out/{wildcards.sample}/STARsolo/Solo.out/GeneFull/soupx {output.H5AD}
            """
        )