#!/bin/bash

PREFETCH_EXEC=$1
BAM2FQ_EXEC=$2
THREADS=$3
SAMPLE=$4
SRR=$(cat ./SRR_list.txt)

cd tmp

# Fetch .bam file
${PREFETCH_EXEC} \
    --verify yes \
    --max-size 999999999999 \
    --force ALL \
    --type all \
    ${SRR}

# --output-directory tmp/ \
# $(cat ./SRR_list.txt)
# --output-file ${SAMPLE}.bam \

# Convert .bam to .fastq's
mkdir -p b2f/
rm -rf b2f/

${BAM2FQ_EXEC} \
    --nthreads ${THREADS} \
    --reads-per-fastq=999999999999999 \
    ${SRR}/*.bam \
    b2f/

# Search for files matching the regex pattern
N_FQs=$(find b2f/ -name *R1*.fastq.gz | wc -l)

# Check if the number of matching files is greater than 0
if [[ ${N_FQs} -gt 1 ]]; then
    #Rename fastqs
    zcat $(find b2f/ -name *R1*.fastq.gz) > merged_R1.fq
    zcat $(find b2f/ -name *R2*.fastq.gz) > merged_R2.fq

    # Compress
    pigz -p${THREADS} merged_*.fq
else
    mv b2f/*/*R1*.fastq.gz merged_R1.fq.gz
    mv b2f/*/*R2*.fastq.gz merged_R2.fq.gz
fi

# Clean up bam2fastq garbage
rm -r b2f/
rm -r ${SRR}/
