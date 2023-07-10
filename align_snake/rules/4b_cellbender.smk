#############################################
## Ambient RNA decontamination
#############################################
# https://github.com/broadinstitute/CellBender/blob/master/docs/source/usage/index.rst
rule ambient_rna_decon_cellbender:
    input:
        # GENEFULL_FILT_MAT = "{DATADIR}/align_out/{sample}/STARsolo/Solo.out/GeneFull/filtered/matrix.mtx.gz",
        # FILT_CELLS = "{DATADIR}/align_out/{sample}/STARsolo/Solo.out/GeneFull/filtered/barcodes.tsv.gz",
        # FILT_FEATS = "{DATADIR}/align_out/{sample}/STARsolo/Solo.out/GeneFull/filtered/features.tsv.gz",
        # GENEFULL_RAW_MAT = "{DATADIR}/align_out/{sample}/STARsolo/Solo.out/GeneFull/raw/matrix.mtx.gz",
        # RAW_CELLS = "{DATADIR}/align_out/{sample}/STARsolo/Solo.out/GeneFull/raw/barcodes.tsv.gz",
        # RAW_FEATS = "{DATADIR}/align_out/{sample}/STARsolo/Solo.out/GeneFull/raw/features.tsv.gz"
        H5AD = DATADIR+"/align_out/{sample}/STARsolo/Solo.out/GeneFull/filtered/decon.h5ad"
    output:
        H5AD = DATADIR+"/align_out/{sample}/STARsolo/Solo.out/GeneFull/filtered/decon.h5ad"
    conda:
        "cellbender"
    threads:
        config["CORES_MID"]
    run:
        script(
            f"""
            cellbender remove-background \
            --input {input.H5AD} \
            --output {output.H5AD} \
            --total-droplets-included 20000 \
            --fpr 0.01 \
            --epochs 150
            """
        )
            # --cuda \
            # --expected-cells 5000 \
