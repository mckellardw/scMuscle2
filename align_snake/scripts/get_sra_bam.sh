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
--output-file ${SRR}.sra \ #TODO
./SRR_list.txt
 # ${SRR}

# Convert .bam to .fastq's
mkdir -p tmp/${SRR}/b2f/
rm -rf tmp/${SRR}/b2f/

${BAM2FQ_EXEC} \
--nthreads ${THREADS} \
--reads-per-fastq=999999999999999 \
${SRA_BAM} \ #TODO
tmp/${SRR}/b2f/

#Rename fastqs
zcat $(find tmp/${SRR}/b2f/ -name *R1*.fastq.gz) > tmp/${SRR}_1.fastq
zcat $(find tmp/${SRR}/b2f/ -name *R2*.fastq.gz) > tmp/${SRR}_2.fastq

# Compress
pigz -p${THREADS} tmp/${SRR}_*.fastq

# Clean up bam2fastq garbage
rm -r tmp/${SRR}/b2f/
