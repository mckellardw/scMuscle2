# `align_snake`
Preprocessing, alignment, QC, and quantification workflow for 10x Genomics single-cell & single-nucleus RNA-seq data (Chromium, v2, v3, or v3.1)
**by David W. McKellar**

#README TODO:
- Required packages & dependencies (add installation via .yml file)
- Write out pipeline details
- Info on sample_sheet format

## **Dependencies:**
- `gget` [v0.2.6](https://github.com/pachterlab/gget)
<!-- - `cutadapt` [v3.4](https://cutadapt.readthedocs.io/en/stable/) -->
- `trim-galore` [v0.6.2]()
- `fastqc` [v##]()
- `STAR` [v2.7.10a](https://github.com/alexdobin/STAR)
- `qualimap` [v.2.2.2a]()
- `pigz` [v2.6](https://zlib.net/pigz/)

*Maybe include*:
- `vsearch` [v#](https://github.com/torognes/vsearch)
- `BLAST`

#### Build `conda` environment
Be sure to [install conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) first!
```
conda create --name align_snake

conda activate align_snake

conda install -c bioconda gget star fastqc trim-galore

conda install  -c conda-forge pigz
```

## **Outputs:**
#TODO: add complete `tree`
