# `align_snake`
Preprocessing, alignment, QC, and quantification workflow for 10x Genomics single-cell & single-nucleus RNA-seq data (Chromium, v2, v3, or v3.1)

#Pipeline TODO:
- Increase data access functionalities
  - Add `aws`/`wget` as a download option
    - Tabula dataset downloads... just `wget` from AWS?
  - Add `local` as a download option under `meta$file.format`
- Add `multiqc` [[link](https://multiqc.info/)] rule to aggregate QC info on alignment/etc.

#README TODO:
- Info on sample_sheet format

## Data to add...
- Tabula Microcebus raw data still not available [link](https://tabula-microcebus.ds.czbiohub.org/whereisthedata)

## **Dependencies:**
- `sra-toolkit` [v3.0.0](https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit)
- `parallel-fastq-dump` [link](https://github.com/rvalieris/parallel-fastq-dump)
- `ffq` [v0.2.1](https://github.com/pachterlab/ffq)
- `bam2fastq` [v1.4.1](https://github.com/10XGenomics/bamtofastq/blob/master/README.md) ([link to download](https://github.com/10XGenomics/bamtofastq/releases))
- `pigz` [v2.6](https://zlib.net/pigz/)
- `pandas` [v##]()
- `itertools` [v##]()
- `gget` [v0.2.6](https://github.com/pachterlab/gget)
- `cutadapt` [v4.1](https://cutadapt.readthedocs.io/en/stable/)
- `fastqc` [v0.11.8](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- `STAR` [v2.7.10b](https://github.com/alexdobin/STAR)
- `qualimap` [v.2.2.2a](http://qualimap.conesalab.org/)
- `multiqc` [v1.13](https://multiqc.info/)


#### **Build `conda` environment:**
Be sure to [install conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) first!
```
conda create --name align_snake

conda activate align_snake

conda install -c bioconda ffq parallel-fastq-dump gget fastqc cutadapt star multiqc kallisto bustools snakemake

conda install itertools
```

## **Pipeiline overview:**
#TODO

## **Runtime info:**
Run snakemake from command line with total desired core usage (`-j num_threads`) and increased latency wait time as well as restarts, because of delays from SRA (`--latency-wait`):
```
conda activate align_snake
snakemake -j 33 --keep-going --max-jobs-per-second 5 --latency 15 --restart-times 3
```
I have found that NCBI does not like when you query their databases too often... The `--max-jobs-per-second` flag should fix this and `--restart-times 3` will retry the download a couple times. If you get a `too many queries` error, just rerun the Snakemake. I also recommend the `--keep-going` or `-k` flag so that the Snakemake doesn't die after encountering the error.


## **Outputs:**
File `tree` for `align_snake` outputs
```
{DATADIR}/{align_out}/{GSM#######}
├── cutadapt_polyA_report.txt
├── cutadapt_polyG_report.txt
├── postTrim_fastqc_R2_out
│   ├── GSM#######_R2_final_fastqc.html
│   └── GSM#######_R2_final_fastqc.zip
├── preTrim_fastqc_R2_out
│   ├── GSM#######_R2_fastqc.html
│   └── GSM#######_R2_fastqc.zip
├── STARsolo
│   ├── Aligned.sortedByCoord.out.bam
│   ├── Log.final.out
│   ├── Log.out
│   ├── Log.progress.out
│   ├── SJ.out.tab
│   ├── Solo.out
│   │   ├── Barcodes.stats
│   │   ├── Gene
│   │   │   ├── Features.stats
│   │   │   ├── filtered
│   │   │   │   ├── barcodes.tsv.gz
│   │   │   │   ├── features.tsv.gz
│   │   │   │   └── matrix.mtx.gz
│   │   │   ├── raw
│   │   │   │   ├── barcodes.tsv.gz
│   │   │   │   ├── features.tsv.gz
│   │   │   │   ├── matrix.mtx.gz
│   │   │   │   └── UniqueAndMult-EM.mtx.gz
│   │   │   ├── Summary.csv
│   │   │   └── UMIperCellSorted.txt
│   │   ├── GeneFull
│   │   │   ├── Features.stats
│   │   │   ├── filtered
│   │   │   │   ├── barcodes.tsv.gz
│   │   │   │   ├── features.tsv.gz
│   │   │   │   └── matrix.mtx.gz
│   │   │   ├── raw
│   │   │   │   ├── barcodes.tsv.gz
│   │   │   │   ├── features.tsv.gz
│   │   │   │   ├── matrix.mtx.gz
│   │   │   │   └── UniqueAndMult-EM.mtx.gz
│   │   │   ├── Summary.csv
│   │   │   └── UMIperCellSorted.txt
│   │   └── Velocyto
│   │       ├── Features.stats
│   │       ├── filtered
│   │       │   ├── ambiguous.mtx.gz
│   │       │   ├── barcodes.tsv.gz
│   │       │   ├── features.tsv.gz
│   │       │   ├── spliced.mtx.gz
│   │       │   └── unspliced.mtx.gz
│   │       ├── raw
│   │       │   ├── ambiguous.mtx.gz
│   │       │   ├── barcodes.tsv.gz
│   │       │   ├── features.tsv.gz
│   │       │   ├── spliced.mtx.gz
│   │       │   └── unspliced.mtx.gz
│   │       └── Summary.csv
│   ├── Unmapped.out.mate1.fastq.gz
│   └── Unmapped.out.mate2.fastq.gz
└── tmp
    ├── GSM#######_R1_final.fq.gz
    └── GSM#######_R2_final.fq.gz
```


## Useful resources on how to upload/download .fastq's from NCBI & other databases:
- ["08. prefetch and fasterq dump", from NCBI](https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump)
- ["Downloading SRA data via GCP", by Tim Stuart](https://timoast.github.io/blog/downloading-sra-data-via-gcp/)
- ["How to use NCBI SRA Toolkit effectively?", by Renesh Bedre](https://www.reneshbedre.com/blog/ncbi_sra_toolkit.html)
