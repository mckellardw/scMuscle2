# **scMuscle2:** single-cell muscle transcriptomics across species
In the first scMuscle project, we attempted to gain a better understanding of the constituent cell types of murine skeletal muscle. In this update, we aim to expand the scope, resolution, and depth of that view, and to improve the availability of our analyses to the greater muscle/single-cell communities.

# How to include your data in scMuscle2
We made a [quick google survey](https://docs.google.com/forms/d/1aueSLPLMHReFEX4Ie1K_o0WEciAznE0HAdOh4eIctVA/prefill) to simplify sharing metadata. Only publicly available data will be included in this study, so be sure to upload to GEO before you do anything else!

Please send any questions to either David McKellar (dwm269@cornell.edu) or to our official scMuscle email address (scmuscle@cornell.edu).

# **TODO:**
- [ ] Finish updating sample list, including new meta data (species, etc)
  - [ ] Add the following metadata columns: age, tissue/muscle, disease/injury status
  - Add HCA data
- [ ] Download fastqs
- [X] Build reference genomes (w/ GENCODE) AND a list of commonly used transgenes
  - [X] Mouse
  - [X] Human
  - [X] ZF
  - [X] rat
  - [X] pig
  - [X] ~~others...~~ Automated this based on the sample metadata! Shout out `gget`
- [ ] Run STARsolo
  - [ ] Add kallisto for comparison?
- [ ] Add ambient RNA decontamination into `align_snake`
  - Tools to try:
    - [ ] SoupX
    - [ ] Cellbender
    - [ ] others...
- [ ] Single-cell analysis (~~Seurat~~ or scanpy?)
- [ ] Benchmark integration w/ [scIB](https://github.com/theislab/scib)? Likely will take a ***LOT*** of compute
- [ ] Figure out multi-species integration (https://github.com/atarashansky/self-assembling-manifold)


# Description of meta data
- **source.label** - shorthand for the publication from which this data was taken (first author's last name plus year final manuscript was published)
- **sample** - sample ID, specific to each lane of a 10x Chromium run (string of characters with no spaces, no periods, and no hyphens/dashes)
- **description** - Brief description of the sample. Not used by any pipeline.
- **include** - Whether or not to include the sample in the analysis (True/False)
- **species** - species from which the sample was derived (Genus species)
- **GSE.accession** - Gene Expression Omnibus (GEO) project accession number (GSE#######)
- **GSM.accession** - GEO sample accession number (GSM#######). *Note* this ID is used for pooling fastq's (SRR #'s) from samples for which authors uploaded multiple files
- **SRA.format** - This column contains the file format in which authors originally uploaded their data. This is used by the `align_snake` pipeline to download and reformat each file so that they may be re-aligned to a common reference. ("fastq"/"bam")
- **SRR.accession** - This column contains the "SRR#########" IDs for each run. These IDs are important because they are used to download the raw sequencing data for each sample (see `align_snake`). For samples which have more than one run, the SRR numbers are delimited with a semi-colon (;). ("SRR#########" or "ERR########")
- **SAMN.accession** - `INSERT_DESCRIPTION_HERE`
- **other.accession** - `INSERT_DESCRIPTION_HERE`
- **source** - This is a citation for each publication
- **chemistry** - All of the samples we collected used the Chromium platform from 10x Genomics. This column specifies what version of the chemistry was used (10x Genomics Chromium v1/v2/v3/v3.1, or other methods). For lots of samples, I used the length of R1 to determine if the sample was v2 or v3 (v3 and v3.1 use the same whitelist, so specifying isn't as important; please let me know if you find an error though!). See `resources/chemistry_sheet.csv` for details on how this is used in alignment (`align_snake`).
- **sequencer** - type of sequencer used
- **sex** - ("M"/"F")
- **age** - measured in... days? weeks? months? #TODO
- **tissue** - #TODO
- **subtissue** - #TODO

## Sources used to find samples
- [NCBI/GEO](https://www.ncbi.nlm.nih.gov/geo/) - keywords used:
  - (muscle) AND (single-cell)
  - (tendon) AND (single-cell)
  - (limb) AND (single-cell)
  - (skeletal) AND (single-cell)
  - (cartilage) AND (single-cell)
- [CNBC/NGDC](https://ngdc.cncb.ac.cn/) - keywords used:
  -
- [10x Genomics - 'Publications'](https://www.10xgenomics.com/resources/publications)
- [panglaoDB](https://panglaodb.se/)
- [Broad Single-Cell Portal](https://singlecell.broadinstitute.org/single_cell)
- [Svensson et al Database](http://www.nxn.se/single-cell-studies/gui)
- [Single-Cell Expression Atlas](https://www.ebi.ac.uk/gxa/sc/home)
- [Human Cell Atlas Data Portal](https://data.humancellatlas.org/explore/projects)

## Other useful tools for exploring sequencing data
- [ffq](https://github.com/pachterlab/ffq) - used to clarify metadata, fill out accession info. Incorporated into the `align_snake` pipeline
- [fastqerq-dump](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump)
- [parallel-fastq-dump](https://github.com/rvalieris/parallel-fastq-dump)

## Single-cell analysis tools
**#TODO**
- [Cellbender](https://github.com/broadinstitute/CellBender) - Ambient RNA removal, [published](https://doi.org/10.1101/791699)
- [Scrublet](https://github.com/swolock/scrublet) - doublet ID and removal, [published](https://www.sciencedirect.com/science/article/pii/S2405471218304745)
- [Scanpy](TODO)
-

# **Workflow**
`align_snake` - snakemake workflow to automate everything between metadata collectino and count matrix preprocessing
1. Download raw sequencing data for all samples (`ffq`,`parallel-fastq-dump`)
2. Build reference genomes (`gget`, `STAR`)
3. Align sequencing data (`STAR`, `kallisto`?)  


`TAR-scRNA-seq`- annotation-free analysis of scRNAseq data [link](https://github.com/fw262/TAR-scRNA-seq)
4. Perform annotation-free analysis `multiTAR` (**#TODO**)


5. Integrate samples for each species separately (`Seurat`/`scanpy`/`Harmony`/etc)
6. Cross-species integrative analysis (`SAMmap`)

## "tissue" annotations
[1] "cell line"                "joint"                   
[3] "limb bud"                 "muscle"                  
[5] "muscle and colon"         "organoid"                
[7] "psc-derived"              "tendon"                  
[9] "tumor"                    ""                        
[11] "#TODO"                    "bone"                    
[13] "cartilage"                "chondrocytes"            
[15] "cranial neural crest"     "ESC"                     
[17] "ESC-derived chondrocytes" "ligament"                
[19] "lung"                     "MSC"                     
[21] "peritoneum"               "skin"                    
[23] "spinal cord"              "xenograft"               
[25] "limb"    
