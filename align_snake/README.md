# `align_snake`
Preprocessing, alignment, QC, and quantification workflow for 10x Genomics single-cell & single-nucleus RNA-seq data (Chromium, v2, v3, or v3.1)
**by David W. McKellar**

#README TODO:
- Write out pipeline details
- Info on sample_sheet format

## **Dependencies:**
- `parallel-fastq-dump` [link](https://github.com/rvalieris/parallel-fastq-dump)
- `ffq` [v0.2.1](https://github.com/pachterlab/ffq)
- `bam2fastq` [v1.4.1](https://github.com/10XGenomics/bamtofastq/blob/master/README.md) ([link to download](https://github.com/10XGenomics/bamtofastq/releases))
- `pigz` [v2.6](https://zlib.net/pigz/)
- `pandas` [v##]()
- `itertools` [v##]()

- `gget` [v0.2.6](https://github.com/pachterlab/gget)
- `cutadapt` [v4.1](https://cutadapt.readthedocs.io/en/stable/)
<!-- - `trim-galore` [v0.6.2](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) -->
- `fastqc` [v0.11.8](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- `STAR` [v2.7.10a](https://github.com/alexdobin/STAR)
- `qualimap` [v.2.2.2a](http://qualimap.conesalab.org/)

*Maybe include*:
- `vsearch` [v#](https://github.com/torognes/vsearch)
- `BLAST`

#### Build `conda` environment
Be sure to [install conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) first!
```
conda create --name align_snake

conda activate align_snake

conda install -c bioconda ffq parallel-fastq-dump gget fastqc cutadapt star

conda install  -c conda-forge pigz pandas

conda install itertools
```

## Pipeiline overview
#TODO

## Runtime info
Run snakemake from command line with total desired core usage (`-j num_threads`) and increased latency wait time as well as restarts, because of delays from SRA (`--latency-wait`):
```
conda activate getfqs_snake
snakemake -j 33 --latency 15 --restart-times 3
```
I have found that NCBI does not like when you query their databases too often, so if you get a `too many queries` error, just rerun the Snakemake.


## **Outputs:**
#TODO: add complete `tree`



## Useful resources on how to upload/download .fastq's from NCBI & other databases:
- ["Downloading SRA data via GCP", by Tim Stuart](https://timoast.github.io/blog/downloading-sra-data-via-gcp/)
- ["How to use NCBI SRA Toolkit effectively?", by Renesh Bedre](https://www.reneshbedre.com/blog/ncbi_sra_toolkit.html)
