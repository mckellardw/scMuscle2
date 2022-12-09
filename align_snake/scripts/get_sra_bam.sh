#!/bin/bash

PREFETCH_EXEC=$1
BAM2FQ_EXEC=$2
THREADS=$3
SRR=$4

# Fetch .bam file
${PREFETCH_EXEC} \
--verify yes \
--max-size 999999999999 \
--force ALL \
--output-directory tmp/ \
--type all \
 ${SRR}

# ./SRR_list.txt
# --output-file ${SRR}.bam \

# Convert .bam to .fastq's
mkdir -p tmp/b2f/
rm -rf tmp/b2f/

${BAM2FQ_EXEC} \
--nthreads ${THREADS} \
--reads-per-fastq=999999999999999 \
${SRR}/*.bam \
tmp/b2f/

# Search for files matching the regex pattern
num_fqs=$(find tmp/b2f/ -name *R1*.fastq.gz | wc -l)

# Check if the number of matching files is greater than 0
if [[ ${num_fqs} -gt 1 ]]; then
    #Rename fastqs
    zcat $(find tmp/b2f/ -name *R1*.fastq.gz) > tmp/merged_1.fastq
    zcat $(find tmp/b2f/ -name *R2*.fastq.gz) > tmp/merged_2.fastq

    # Compress
    pigz -p${THREADS} tmp/${SRR}_*.fastq
else
    mv tmp/b2f/*/*R1*.fastq.gz tmp/merged_1.fastq.gz
    mv tmp/b2f/*/*R2*.fastq.gz tmp/merged_2.fastq.gz
fi

# Clean up bam2fastq garbage
rm -r tmp/b2f/
rm -r tmp/${SRR}/
