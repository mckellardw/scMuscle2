# Doublet & UMI count filtering
rule qc_preprocessing:
    input:
        H5AD = DATADIR+"/align_out/{sample}/STARsolo/Solo.out/GeneFull/soupx/decon.h5ad"
    output:
        H5AD = DATADIR+"/align_out/{sample}/STARsolo/Solo.out/GeneFull/soupx/deconQC.h5ad"
    params:
        min_genes = 500,
        min_counts = 1000,
        max_mito = 40
    log:
        DATADIR+"/align_out/{sample}/preprocessing.log"
    threads:
        1
    run:
        shell(
            f"""
            python scripts/qc_h5ad.py {input.H5AD} {output.H5AD} {params.min_genes} {params.min_counts} {params.max_mito} 2> {log}
            """
        )
