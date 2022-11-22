# initialize & cache anndata object
rule cache_h5ad:
    input:
        GENEFULLMAT = DATADIR+"/align_out/{sample}/STARsolo/Solo.out/GeneFull/filtered/matrix.mtx.gz"
    output:
        H5AD = "{CACHEDIR}/{sample}.h5ad"
    threads:
        1
        # config["CORES_LO"]
    run:
        # from scanpy import read_10x_mtx#, settings
        import scanpy
        # settings.cachedir = CACHEDIR
        adata = scanpy.read_10x_mtx(
            path="{DATADIR}/align_out/{sample}/STARsolo/Solo.out/GeneFull/filtered",
            var_names='gene_symbols',
            make_unique=True,
            cache=False
        )

        adata.write(output.H5AD)
