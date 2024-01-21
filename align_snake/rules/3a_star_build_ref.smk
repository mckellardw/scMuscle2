# Build reference genomes for STARsolo alignment/quantification
#TODO: gget import error for some reason?
rule build_refs:
    output:
        # REF_METADATA = expand("{REFDIR}/{SPECIES}/STAR/metadata.json", REFDIR=config["REFDIR"], SPECIES=SPECIES),
        REF = expand("{REFDIR}/{SPECIES}/STAR/Genome", REFDIR=REFDIR, SPECIES=SPECIES) # Reference genomes
    threads:
        config["CORES_HI"]
    log:
        f"{REFDIR}/star.log"
    run:
        #TODO- import error from gget python module, although the command line tool works? Seems like it is a
        # from gget import ref
        # available_species = ref(species="NA", list_species=True)
        shell(
            f"{EXEC['GGET']} ref -l > resources/gget_species.txt"
        )

        from pandas import read_csv
        available_species = read_csv("resources/gget_species.txt",header=None)[0].values.tolist()

        for S in SPECIES:
            if S in available_species:
                print(f"{REFDIR}/{S}/STAR/Genome")
                if not path.exists(f"{REFDIR}/{S}/STAR/Genome"):
                    print(f"Downloading genome sequence and annotations for {S} to {REFDIR}/{S}")
                    shell(
                        f"""
                        mkdir -p {REFDIR}/{S}
                        cd {REFDIR}/{S}

                        {EXEC['GGET']} ref \
                        --out {REFDIR}/{S}/metadata.json \
                        --which gtf,dna \
                        --download \
                        {S}

                        gunzip {REFDIR}/{S}/*.gz
                        """
                    )

                    # Build reference for STAR
                    print(f"Building STAR reference for {S}...\n")
                    shell(
                        f"""
                        {EXEC['STAR']} \
                        --runThreadN {threads} \
                        --runMode genomeGenerate \
                        --genomeDir {REFDIR}/{S}/STAR \
                        --genomeFastaFiles $(ls -t {REFDIR}/{S}/*.fa) \
                        --sjdbGTFfile $(ls -t {REFDIR}/{S}/*.gtf) \
                        --sjdbGTFfeatureExon exon

                        pigz -p {threads} $(ls -t {REFDIR}/{S}/*.fa)
                        pigz -p {threads} $(ls -t {REFDIR}/{S}/*.gtf)
                        """
                    )
                else:
                    print(f"STAR reference for {S} already found...")
            else:
                print(f"ENSEMBL reference for {S} not found by gget... :(")
