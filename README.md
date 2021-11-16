# **scMuscle2:** single-cell muscle transcriptomics across species
In the scMuscle project, we attempted to gain a better understanding of the constituent cell types of murine skeletal muscle. In this update, we aim to expand the scope and resolution of that view, and to improve the availability of our analyses to the greater muscle/single-cell communities.

## TODO:

- Finish updating sample list, including new meta data (species, etc) [ ]
- Download fastqs [ ]
- Build reference genomes (w/ GENCODE)
  - Mouse [X]
  - Human [X]
  - ZF [ ]
- Run STARsolo [ ]
  - Add species selection step to STARsolo pipeline
- Single-cell analysis (seurat or scanpy?) [ ]
- Figure out multi-species integration (https://github.com/atarashansky/self-assembling-manifold) [ ]


## Description of meta data
- **source.label** - shorthand for the publication from which this data was taken (first author's last name plus year final manuscript was published)
- **sample** - sample ID, specific to each lane of a 10x Chromium run
- **description** - INSERT_DESCRIPTION_HERE
- **include** - INSERT_DESCRIPTION_HERE
- **species** - species from which the sample was derived (Full genus and species) ***#TODO update this in metadata***
- **GSE.accession** - INSERT_DESCRIPTION_HERE
- **GSM.accession** - INSERT_DESCRIPTION_HERE
- **SRA.is.fastq** - INSERT_DESCRIPTION_HERE
- **SRR.accession** - INSERT_DESCRIPTION_HERE
- **SAMN.accession** - INSERT_DESCRIPTION_HERE
- **other.accession** - INSERT_DESCRIPTION_HERE
- **source** - INSERT_DESCRIPTION_HERE
-

## Sources used to find samples
- [NCBI/GEO](https://www.ncbi.nlm.nih.gov/geo/) - keywords used: ""
- [10x Genomics - Publications](https://www.10xgenomics.com/resources/publications) -
- [panglaoDB](https://panglaodb.se/)

## Other useful tools for exploring sequencing data
- [ffq](https://github.com/pachterlab/ffq)
- [fastqerq-dump](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump)
- [parallel-fastq-dump]()
