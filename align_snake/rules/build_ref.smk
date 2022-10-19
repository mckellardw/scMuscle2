# Build reference genomes for STARsolo alignment/quantification
#TODO: gget import error for some reason?
rule build_refs:
    output:
        # REF_METADATA = expand("{REFDIR}/{SPECIES}/STAR/metadata.json", REFDIR=config["REFDIR"], SPECIES=SPECIES),
        REF = expand("{REFDIR}/{SPECIES}/STAR/SAindex", REFDIR=config["REFDIR"], SPECIES=SPECIES) # Reference genomes
    threads:
        config["CORES_HI"]
    run:
        from gget import ref

        available_species = ref(species="NA", list_species=True)
        for S in SPECIES:
            if S in available_species:
                if path.isfile(f"{REFDIR}/{S}/STAR/SAindex"):
                    print(f"Downloading genome sequence and annotations for {S} to {REFDIR}/{S}")

                    shell(
                        f"""
                        mkdir -p {REFDIR}/{wildcards.SPECIES}
                        cd {REFDIR}/{wildcards.SPECIES}

                        {GGET_EXEC} ref \
                        --out {REFDIR}/{wildcards.SPECIES}/metadata.json \
                        --which gtf,dna \
                        --download \
                        {wildcards.SPECIES}

                        gunzip {REFDIR}/{wildcards.SPECIES}/*.gz
                        """
                    )

                    # Build reference for STAR
                    print(f"Building STAR reference for {wildcards.SPECIES}...\n")
                    shell(
                        f"""
                        {STAR_EXEC} \
                        --runThreadN {threads} \
                        --runMode genomeGenerate \
                        --genomeDir {REFDIR}/{wildcards.SPECIES}/STAR \
                        --genomeFastaFiles $(ls -t {REFDIR}/{wildcards.SPECIES}/*.fa) \
                        --sjdbGTFfile $(ls -t {REFDIR}/{wildcards.SPECIES}/*.gtf) \
                        --sjdbGTFfeatureExon exon

                        pigz -p {threads} $(ls -t {REFDIR}/{wildcards.SPECIES}/*.fa)
                        pigz -p {threads} $(ls -t {REFDIR}/{wildcards.SPECIES}/*.gtf)
                        """
                    )
                else:
                    print(f"STAR reference for {S} already found...")
            else:
                print(f"ENSEMBL reference for {S} not found by gget... :(")
