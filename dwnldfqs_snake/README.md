# download_fastqs_snakemake

## Runtime info
Run snakemake from command line with total desired core usage (`-j num_threads`) and increased latency wait time as well as restarts, because of delays from SRA (`--latency-wait`):

```
snakemake -j 33 --latency 15 --restart-times 3
```

I have also found that NCBI does not like when you query their databases too often, so if you get a `too many queries` error, just rerun the workflow.

## Requirements
- `parallel-fastq-dump` [link](https://github.com/rvalieris/parallel-fastq-dump)
- `ffq` [link](https://github.com/pachterlab/ffq)
- `bam2fastq` [link](https://github.com/10XGenomics/bamtofastq/blob/master/README.md) (built w/ v1.4.1)
- `pigz` (parallelized gzip, huge time saver)
- `pandas`
- `itertools`


### Useful resources on how to upload/download .fastq's from NCBI & other databases:
- ["Downloading SRA data via GCP", by Tim Stuart](https://timoast.github.io/blog/downloading-sra-data-via-gcp/)
- ["How to use NCBI SRA Toolkit effectively?", by Renesh Bedre](https://www.reneshbedre.com/blog/ncbi_sra_toolkit.html)
