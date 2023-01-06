# initialize & cache anndata object
rule cache_preQC_h5ad:
    input:
        GENEFULLMAT = DATADIR+"/align_out/{sample}/STARsolo/Solo.out/GeneFull/filtered/matrix.mtx.gz"
    output:
        H5AD = DATADIR+"/align_out/{sample}/STARsolo/Solo.out/GeneFull/filtered/preDecon.h5ad"
    threads:
        1
    run:
        shell(
            f"""
            python scripts/cache_h5ad.py {DATADIR}/align_out/{wildcards.sample}/STARsolo/Solo.out/GeneFull/filtered {output.H5AD}
            """
        )