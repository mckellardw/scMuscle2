
# Description of metadata
#TODO- update to match cellxgene requirements

#### Source/citation info:
- **source.label** - shorthand for the publication from which this data was taken (first author's last name plus year final manuscript was published)
- **sample** - sample ID, specific to each lane of a 10x Chromium run (string of characters with no spaces, no periods, and no hyphens/dashes)
- **description** - Brief description of the sample. Not used by any pipeline.
- **include** - Whether or not to include the sample in the analysis (True/False)
- **source** - This is an abbreviate citation for each publication (Author et al, Journal, year)

#### Biological metadata & experimental conditions
- **species** - species from which the sample was derived (Genus species)
- **sex** - ("M"/"F")
- **age** - measured in... days? weeks? months? #TODO
- **tissue** - #TODO
- **subtissue** - #TODO

#### Technical info:
- **chemistry** - This column specifies what chemistry/method was used to prepare the data (10x Genomics Chromium v1/v2/v3/v3.1, or other methods). ***Note*** For lots of samples where chemistry was not stated beyond "10x Genomics", I used the length of R1 to determine if the sample was v2 or v3 (v3 and v3.1 use the same R1 length, so specifying isn't as important; please let me know if you find an error though!). See `resources/chemistry_sheet.csv` for details on how this is used in alignment (`align_snake`).
- **sequencer** - type of sequencer used

#### Accession info:
- **file.format** - This column contains the file format in which authors originally uploaded their data. This is used by the `align_snake` pipeline to download and reformat each file so that they may be re-aligned to a common reference. Options:
  - "fastq"
  - "bam"
  - "aws" - file is downloaded w/ `wget` and renamed to match the accession info `#TODO`
  - "local" - option for files that you have downloaded manually, or files that aren't public (yet...) `#TODO`
- **GSE.accession** - Gene Expression Omnibus (GEO) project accession number (GSE#######)
- **GSM.accession** - GEO sample accession number (GSM#######). *Note* this ID is used for pooling fastq's (SRR #'s) from samples for which authors uploaded multiple files
- **SRR.accession** - This column contains the "SRR#########" IDs for each run. These IDs are important because they are used to download the raw sequencing data for each sample (see `align_snake`). For samples which have more than one run, the SRR numbers are delimited with a semi-colon (;). ("SRR#########" or "ERR########")
- **SAMN.accession** - `INSERT_DESCRIPTION_HERE`
- **other.accession** - `INSERT_DESCRIPTION_HERE`
- **file.link** - link that can be used to directly download files with `wget`
