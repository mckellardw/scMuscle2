# Doublet & UMI count filtering
rule cache_h5ad:
    input:
        H5AD = DATADIR+"/align_out/{sample}/STARsolo/Solo.out/GeneFull/filtered/decon.h5ad"
    output:
        H5AD = DATADIR+"/align_out/{sample}/STARsolo/Solo.out/GeneFull/filtered/deconQC.h5ad"
    params:
        min_genes = 500,
        min_counts = 1000,
        max_mito = 40
    threads:
        1
    run:
        shell(
            f"""
            mkdir -p {CACHEDIR}/qc/
            python scripts/qc_h5ad.py {input.H5AD} {output.H5AD} {params.min_genes} {params.min_counts} {params.max_mito}
            """
        )
