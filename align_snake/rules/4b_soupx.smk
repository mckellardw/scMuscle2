#############################################
## Ambient RNA decontamination
#############################################
# Run soupx on the STARsolo outputs. 
## *Note-* the features file for the soupx matrix is copied over from the `.../filtered/` directory to avoid anndata/scanpy errors. They are identical, except for the  number of columns
rule ambient_rna_decon_soupx:
    input:
        GENEFULL_FILT_MAT = DATADIR+"/align_out/{sample}/STARsolo/Solo.out/GeneFull/filtered/matrix.mtx.gz",
        GENEFULL_FILT_CELLS = DATADIR+"/align_out/{sample}/STARsolo/Solo.out/GeneFull/filtered/barcodes.tsv.gz",
        GENEFULL_FILT_FEATS = DATADIR+"/align_out/{sample}/STARsolo/Solo.out/GeneFull/filtered/features.tsv.gz",
        GENEFULL_RAW_MAT = DATADIR+"/align_out/{sample}/STARsolo/Solo.out/GeneFull/raw/matrix.mtx.gz",
        RAW_CELLS = DATADIR+"/align_out/{sample}/STARsolo/Solo.out/GeneFull/raw/barcodes.tsv.gz",
        RAW_FEATS = DATADIR+"/align_out/{sample}/STARsolo/Solo.out/GeneFull/raw/features.tsv.gz"
    output:
        SOUPX_MAT = DATADIR+"/align_out/{sample}/STARsolo/Solo.out/GeneFull/soupx/matrix.mtx.gz",
        SOUPX_CELLS = DATADIR+"/align_out/{sample}/STARsolo/Solo.out/GeneFull/soupx/barcodes.tsv.gz",
        SOUPX_FEATS = DATADIR+"/align_out/{sample}/STARsolo/Solo.out/GeneFull/soupx/features.tsv.gz"
    log:
        DATADIR+"/align_out/{sample}/soupx.log"
    threads:
        1
        # config["CORES_LO"] 
    run:
        run_decon = META.loc[META["GSM.accession"] == wildcards.sample, "ambient.decon"].values[0]
        if run_decon:
            shell(
                f"""
                cd {DATADIR}/align_out/{wildcards.sample}
                Rscript {PRODIR}/align_snake/scripts/starsolo_soupx.R {DATADIR}/align_out/{wildcards.sample}/STARsolo/Solo.out/GeneFull {threads} 2> {log}
                cp {input.GENEFULL_FILT_FEATS} {output.SOUPX_FEATS}
                """
            )
        else:
            shell(
                f"""                
                mkdir -p {DATADIR}/align_out/{wildcards.sample}/STARsolo/Solo.out/GeneFull/soupx
                cp {input.GENEFULL_FILT_MAT} {output.SOUPX_MAT}
                cp {input.GENEFULL_FILT_CELLS} {output.SOUPX_CELLS}
                cp {input.GENEFULL_FILT_FEATS} {output.SOUPX_FEATS}

                echo "Copied GeneFull matrix over to soupx folder, no ambient RNA decontamination run." > {log}
                """
            )
            
        shell(
            f"ls -R {DATADIR}/align_out/{wildcards.sample}/STARsolo/Solo.out/GeneFull/ >> {log}"
        )

# initialize & cache the **soupx** counts as an anndata file for easier loading later
rule cache_soupx_h5ad:
    input:
        GENEFULLMAT = DATADIR+"/align_out/{sample}/STARsolo/Solo.out/GeneFull/soupx/matrix.mtx.gz"
    output:
        H5AD = DATADIR+"/align_out/{sample}/STARsolo/Solo.out/GeneFull/soupx/decon.h5ad"
    params:
        var_names = "gene_symbols" # scanpy.read_10x_mtx(), index for features is zero (1st column)
    threads:
        1
    run:
        shell(
            f"""
            python {PRODIR}/align_snake/scripts/cache_h5ad.py {DATADIR}/align_out/{wildcards.sample}/STARsolo/Solo.out/GeneFull/soupx {output.H5AD} {params.var_names}
            """
        )