#############################################
## Get metadata with `ffq`
#############################################
# Download metadata for each individual SRR (run) ID
## Info: https://github.com/pachterlab/ffq
rule get_metadata_ffq:
    output:
        METAJSON = "{METADIR}/ffq/{SRR}.json"
    threads:
        config["CORES_MID"] #Limit number of concurrent 
    log:
        "{METADIR}/logs/{SRR}.log"
    run:
        if wildcards.SRR != "NA":
            shell(
                f"""
                {EXEC['FFQ']} -o {output.METAJSON} {wildcards.SRR} \
                2> {log}
                """
            )
        else:
            shell( # For samples which are missing an SRR ID...
                f"""
                touch {output.METAJSON}
                """
            )

# Combine ffq-downloaded .json"s into a big csv for parsing
#TODO- fix this...
rule merge_metadata_ffq:
    input:
        expand("{METADIR}/ffq/{SRR}.json", METADIR=METADIR, SRR=SRR_LIST)
    output:
        METACSV = "{METADIR}/merged_metadata_ffq.csv"
    run:
        df_list = list()
        # merged_df = pd.DataFrame()
        for f in input: # load metadata files for each SRR
            if path.basename(f) == "NA.json":
                pass # ignore NAs
            else:
                with open(f) as data_file:
                    data = json.load(data_file)
                df = pd.json_normalize(data).T
                tmp = df.index[0].split(".")[0] + "." # Get SRR ID...
                df = df.rename(index = lambda x: x.strip(tmp)) # Remove SRR ID from row names
                df_list.append(df.T) # transpose so that each row is a different run

        print("Concatenating metadata...")
        df_out = pd.concat(df_list)

        df_out.to_csv(
            output.METACSV, 
            index=False
        )

#-------------------------------------

rule get_metadata_fastqdl:
    output:
        META = "{METADIR}/fastqdl/{SRR}-run-info.tsv"
    params:
        SLEEP = 14,
        MAX_ATTEMPTS = 4
    threads:
        config["CORES_LO"] #Limit number of concurrent 
    log:
        "{METADIR}/logs/{SRR}_fastqdl.log"
    run:
        if wildcards.SRR != "NA":
            shell(
                f"""
                {EXEC['FASTQDL']} \
                --accession {wildcards.SRR} \
                --sleep {params.SLEEP} \
                --max-attempts {params.MAX_ATTEMPTS} \
                --only-download-metadata \
                --prefix {wildcards.SRR} \
                --outdir $(dirname {output.META}) \
                2> {log}
                """
            )
        else:
            shell( # For samples which are missing an SRR ID...
                f"""
                touch {output.META}
                """
            )

rule merge_metadata_fastqdl:
    input:
        expand("{METADIR}/fastqdl/{SRR}-run-info.tsv", METADIR=METADIR, SRR=SRR_LIST)
    output:
        METACSV = "{METADIR}/merged_metadata_fastqdl.csv"
    run:
        df_list = list()
        # merged_df = pd.DataFrame()
        for f in input: # load metadata files for each SRR
            if path.basename(f) == "NA.json":
                pass # ignore NAs
            else:
                with open(f) as data_file:
                    data = json.load(data_file)
                df = pd.json_normalize(data).T
                tmp = df.index[0].split(".")[0] + "." # Get SRR ID...
                df = df.rename(index = lambda x: x.strip(tmp)) # Remove SRR ID from row names
                df_list.append(df.T) # transpose so that each row is a different run

        print("Concatenating metadata...")
        df_out = pd.concat(df_list)

        df_out.to_csv(
            output.METACSV, 
            index=False
        )