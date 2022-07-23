# `dwnldfqs_snake`

## Runtime info
Run snakemake from command line with total desired core usage (`-j num_threads`) and increased latency wait time as well as restarts, because of delays from SRA (`--latency-wait`):

`snakemake -j 33 --latency 15 --restart-times 3`

I have found that NCBI does not like when you query their databases too often, so if you get a `too many queries` error, just rerun the Snakemake.

## **Dependencies:**
- `parallel-fastq-dump` [link](https://github.com/rvalieris/parallel-fastq-dump)
- `ffq` [v0.2.1](https://github.com/pachterlab/ffq)
- `bam2fastq` [v1.4.1](https://github.com/10XGenomics/bamtofastq/blob/master/README.md) ([link to download](https://github.com/10XGenomics/bamtofastq/releases))
- `pigz`[v2.6]()
- `pandas` [v##]()
- `itertools` [v##]()

#### Build `conda` environment
Be sure to [install conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) first!
```
conda create --name getfqs_snake

conda activate getfqs_snake

conda install -c bioconda gget ffq parallel-fastq-dump

conda install -c conda-forge pigz pandas itertools
```

## **Outputs:**
#TODO

### Useful resources on how to upload/download .fastq's from NCBI & other databases:
- ["Downloading SRA data via GCP", by Tim Stuart](https://timoast.github.io/blog/downloading-sra-data-via-gcp/)
- ["How to use NCBI SRA Toolkit effectively?", by Renesh Bedre](https://www.reneshbedre.com/blog/ncbi_sra_toolkit.html)
