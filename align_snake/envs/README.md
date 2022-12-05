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
