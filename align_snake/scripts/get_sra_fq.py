import sys
# import path
from os import path, system, rename
import glob

#TODO- change this from positional arguments...
args = sys.argv
PREFETCH_EXEC = args[1]
FQD_EXEC = args[2]
THREADS = args[3]
MEMLIMIT = args[4]

# Read in SRR_list
SRR_file = open("SRR_list.txt", "r")
data = SRR_file.read()
SRR_list = data.split("\n")
SRR_file.close()

SRR_list = [i for i in SRR_list if i] #remove empty lines
SRR_list = [i.replace(" ","") for i in SRR_list]# remove space characters

for SRR in SRR_list:
    system(
        f"""
        {PREFETCH_EXEC} \
        --verify yes \
        --type sra \
        --max-size 999G \
        --output-file tmp/{SRR}.sra \
        {SRR}

        {FQD_EXEC} \
        --threads {THREADS} \
        --mem {MEMLIMIT} \
        --outdir tmp/ \
        --temp tmp/ \
        --split-files \
        --include-technical \
        tmp/{SRR}.sra && \
        rm tmp/{SRR}.sra
        """
    )

    # Remove index reads (readlength=8 or readlength=10) | readlength < 26 (just using 20bp below...)
    system(
        f"""
        for SRR_FQ in tmp/{SRR}*.fastq
        do
            if [[ $(cat $SRR_FQ | awk 'NR==2 {{ print }}' | wc -m) -lt 27 ]]; then
                rm $SRR_FQ
            fi
        done
        """
    )

    # Rename fastqs so that R1 and R2 files can be merged
    fq_files = glob.glob(f"tmp/{SRR}*.fastq")

    sra_sizes = {}
    for fq_file in fq_files:
        sra_sizes[fq_file] = path.getsize(fq_file)

    rename_dict = {
        "R2": max(sra_sizes),
        "R1": min(sra_sizes)
    }

    rename(rename_dict["R1"], f"tmp/{SRR}_R1.fastq")
    rename(rename_dict["R2"], f"tmp/{SRR}_R2.fastq")
#end of loop

# merge .fastqs into GSM####.fastq, then remove individual .fastq's
system(
    f"""
    cat tmp/*_R1.fastq > tmp/merged_R1.fq
    cat tmp/*_R2.fastq > tmp/merged_R2.fq

    rm tmp/SRR*.fastq
    """
)

# compress
system(
    f"""
    pigz -p{THREADS} --force tmp/merged_R*.fq
    """
)
