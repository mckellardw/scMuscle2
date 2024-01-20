# Aggregate species-level atlases
rule atlas_init:
    input:
        expand( 
            "{DATADIR}/align_out/{sample}/STARsolo/Solo.out/GeneFull/soupx/deconQC_meta.h5ad", 
            DATADIR=DATADIR,
            sample=SAMPLES # GSM.accession
        ), 
    output:
        H5AD = expand(
            "{DATADIR}/atlas/{SPECIES}.h5ad", 
            DATADIR=DATADIR,
            SPECIES=SPECIES
        )
    params:
        min_genes = 500,
        min_counts = 1000,
        max_mito = 40
    # log:
    #     DATADIR+"/atlas/{species}.merging.log"
    threads:
        1
    run:
        for S in SPECIES:
            tmp_samples = {gsm_id for gsm_id, species in SPECIES_DICT.items() if species == S}
            ad_list = [f"{DATADIR}/align_out/{tmp_sample}/STARsolo/Solo.out/GeneFull/soupx/deconQC_meta.h5ad" for tmp_sample in tmp_samples]
            shell(
                f"""
                python scripts/merge_h5ads.py \
                    "{DATADIR}/atlas/{S}.h5ad" \
                    "GSM.accession" \
                    {" ".join(ad_list)} \
                    | tee {DATADIR+"/atlas/{S}.merging.log"}
                """
            )
