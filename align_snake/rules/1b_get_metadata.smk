#############################################
## Get metadata with `ffq`
#############################################
# Download metadata for each individual SRR (run) ID
## Info: https://github.com/pachterlab/ffq
rule get_metadata:
    output:
        METAJSON = "{METADIR}/{SRR}.json"
    threads:
        config["CORES_MID"] #Limit number of concurrent 
    run:
        if wildcards.SRR != "NA":
            shell(
                f"""
                {EXEC['FFQ']} -o {output.METAJSON} {wildcards.SRR}
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
rule merge_metadata:
    input:
        expand("{METADIR}/{SRR}.json", METADIR=METADIR, SRR=SRR_LIST)
    output:
        METACSV = "{PRODIR}/merged_metadata.csv"
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
