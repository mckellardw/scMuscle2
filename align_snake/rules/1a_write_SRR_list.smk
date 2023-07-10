# Save the list of SRR ids for each sample as a .txt file
rule write_SRR_list:
    output:
        SRR_LIST = "{DATADIR}/align_out/{sample}/SRR_list.txt"
    run:
        SRR_list = SRR_DICT[wildcards.sample]
        with open(output.SRR_LIST, "w") as f:
            f.write('\n'.join(SRR_list))
