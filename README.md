# **scMuscle2:** single-cell muscle transcriptomics across species
In the scMuscle project, we attempted to gain a better understanding of the constituent cell types of murine skeletal muscle. In this update, we aim to expand the scope and resolution of that view, and to improve the availability of our analyses to the greater muscle/single-cell communities.

# TODO:
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


# Description of meta data
- **source.label** - shorthand for the publication from which this data was taken (first author's last name plus year final manuscript was published)
- **sample** - sample ID, specific to each lane of a 10x Chromium run (string of characters with no spaces, no periods, and no hyphens/dashes)
- **description** - Brief description of the sample. Not used by any pipeline.
- **include** - Whether or not to include the sample in the analysis (T/F)
- **species** - species from which the sample was derived (Genus species)
- **GSE.accession** - INSERT_DESCRIPTION_HERE
- **GSM.accession** - INSERT_DESCRIPTION_HERE
- **SRA.format** - This column contains the file format in which authors originally uploaded their data. This is used by the `dwnldqs_snake` pipeline to download and reformat each file so that they may be re-aligned to a common reference. (fastq/bam)
- **SRR.accession** - This column contains the "SRR#########" IDs for each run. These IDs are important because they are used to download the raw sequencing data for each sample (see `dwnldfqs_snake`). For samples which have more than one run, the SRR numbers are delimited with a semi-colon. (SRR#########)
- **SAMN.accession** - INSERT_DESCRIPTION_HERE
- **wget_link** - Link to download the original file, simply using `wget`. *Important to note that this is needed for any data for which the reads were not originally uploaded as .fastq files*
- **other.accession** - INSERT_DESCRIPTION_HERE
- **source** - This is a citation for each publication
- **chemistry** - All of the samples we collected used the Chromium platform from 10x Genomics. This column specifies what version of the chemistry was used. (v2/v3/v3.1)

## Sources used to find samples
- [NCBI/GEO](https://www.ncbi.nlm.nih.gov/geo/) - keywords used: ""
- [10x Genomics - Publications](https://www.10xgenomics.com/resources/publications) -
- [panglaoDB](https://panglaodb.se/)

## Other useful tools for exploring sequencing data
- [ffq](https://github.com/pachterlab/ffq)
- [fastqerq-dump](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump)
- [parallel-fastq-dump](https://github.com/rvalieris/parallel-fastq-dump)

# Workflow

1. Build reference genomes `bldrefs_snake`
2. Download raw sequencing data for all samples `dwnldfqs_snake`
3. Align sequencing data `solo_snake`
4. Perform annotation-free analysis `multiTAR`
5. Integrate samples for each species separately (`Seurat`/`scanpy`/`Harmony`/etc)
6. Cross-species integrative analysis (`SAMmap`)
