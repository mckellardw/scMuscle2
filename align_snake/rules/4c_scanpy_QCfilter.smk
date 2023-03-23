# Doublet & UMI count filtering
rule qc_preprocessing:
    input:
        H5AD = DATADIR+"/align_out/{sample}/STARsolo/Solo.out/GeneFull/soupx/decon.h5ad"
    output:
        H5AD = temp(DATADIR+"/align_out/{sample}/STARsolo/Solo.out/GeneFull/soupx/deconQC.h5ad")
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


# Add metadata to the h5ad's
rule add_meta_to_h5ad:
    input:
        H5AD = DATADIR+"/align_out/{sample}/STARsolo/Solo.out/GeneFull/soupx/deconQC.h5ad"
    output:
        H5AD = DATADIR+"/align_out/{sample}/STARsolo/Solo.out/GeneFull/soupx/deconQC_meta.h5ad"
    params:
    threads:
        1
    run:
        # Subset metadata for just this sample
        tmp_meta = META[list(META["GSM.accession"]==wildcards.sample)]

        # Write the metadata as a .csv in the sample directory
        tmp_meta.to_csv(f"{DATADIR}/align_out/{wildcards.sample}/meta.csv")

        # Run script to add metadata.csv to the h5ad, and save the updated object
        shell(
            f"""
            python scripts/add_meta_h5ad.py {input.H5AD} {output.H5AD} {DATADIR}/align_out/{wildcards.sample}/meta.csv
            """
        )
