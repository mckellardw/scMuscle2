# download_fastqs_snakemake

## Runtime info
Run snakemake from command line with total desired core usage (`-j num_threads`) and increased latency wait time as well as restarts, because of delays from SRA (`--latency-wait`):

`snakemake -j 33 --latency 15 --restart-times 3`


## Requirements
- `parallel-fastq-dump` [link](https://github.com/rvalieris/parallel-fastq-dump)
- `pandas`
- `itertools`

- `ffq` [link](TODO)
- ``
