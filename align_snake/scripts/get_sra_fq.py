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
print(SRR_list)

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

    # Remove technical reads (readlength<10)
    system(
        f"""
        for SRR_FQ in tmp/{SRR}*.fastq
        do
            if cat $SRR_FQ | head | grep -q "length=8"; then
                rm $SRR_FQ
            fi
        done
        """
    )

    # Rename fastqs
    fq_files = glob.glob(f"tmp/{SRR}*.fastq")

    sra_sizes = {}
    for fq_file in fq_files:
        sra_sizes[fq_file] = path.getsize(fq_file)

    print(sra_sizes)
    rename_dict = {
        "R2": max(sra_sizes),
        "R1": min(sra_sizes)
    }

    rename(rename_dict["R1"], f"tmp/{SRR}_R1.fastq")
    rename(rename_dict["R2"], f"tmp/{SRR}_R2.fastq")
#end of loop

# merge fastqs into GSM####.fastq, then remove individual fastq's
system(
    f"""
    cat tmp/*_R1.fastq > tmp/merged_R1.fastq
    cat tmp/*_R2.fastq > tmp/merged_R2.fastq

    rm tmp/SRR*.fastq
    """
)

# compress
system(
    f"""
    pigz -p{THREADS} --force tmp/merged_R*.fastq
    """
)