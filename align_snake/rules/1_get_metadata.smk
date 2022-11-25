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
    shell:
        """
        {FFQ_EXEC} -o {output.METAJSON} {wildcards.SRR}
        """

# Combine ffq-downloaded .json"s into a big csv for parsing
#TODO- fix this...
rule merge_metadata:
    input:
        expand("{METADIR}/{SRR}.json", METADIR=config["METADIR"], SRR=SRR_LIST)
    output:
        METACSV = "{METADIR}/merged_metadata.csv"
    run:
        df_list = list()
        for f in input: # load metadata files for each SRR
            with open(f) as data_file:
                data = json.load(data_file)
            df = pd.json_normalize(data).T
            tmp = df.index[0].split(".")[0] + "." # Get SRR ID...
            df = df.rename(index = lambda x: x.strip(tmp)) # Remove SRR ID from row names
            df_list.append(df.T) # transpose so that each row is a different run

        print("Concatenating metadata...")
        df_out = pd.concat(df_list)

        df_out.to_csv(output.METACSV, index=False)
