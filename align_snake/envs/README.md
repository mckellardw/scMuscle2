# envs
conda environments used in `align_snake`

Two copies of each environment are included, one with the specific build we used (`*-builds.yml`) and the other which has the package versions rather than the specific build.

Example code for writing the .yml files:
```
# cd /path/to/scMuscle2/align_snake
conda env export > envs/align_snake-builds.yml
conda env export --no-builds > envs/align_snake-versions.yml
```

To build the conda environment from scratch, try the following command:
```
conda env create --file envs/align_snake-versions.yml
```

Or to do it a bit quicker...
```
conda  create --name scm2 -c bioconda star==2.7.10b cutadapt==4.1 gget ffq 

conda activate scm2
conda install -c conda-forge pigz pandas itertools scanpy
conda install -c bioconda fastqc multiqc qualimap samtools
```