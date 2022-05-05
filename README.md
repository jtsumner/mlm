# microbiome-snakemake
This is a microbiome snakemake workflow to identify and analyze viral contigs from the built miroenvironment

Note: This repository was built using https://github.com/snakemake-workflows/cookiecutter-snakemake-workflow as a temlate

-------------

## Notes on snakemake 

This pipeline has been tested using snakemake and mamba installed in a single conda environment. Software versions:
* snakemake version 7.3.8 
* mamba 0.15.3
* conda 4.10.3


Execute to create a DAG visualization of the pipeline

```
snakemake --forceall --dag | dot -Tpdf > dag.pdf
```

Execute to create a rule graph visualization 

```
snakemake --forceall --rulegraph | dot -Tpdf > dag.pdf
```


calculate_unifrac.R from 
https://github.com/biobakery/MetaPhlAn/blob/master/metaphlan/utils/calculate_unifrac.R

-------------

## Installation

1. Create a new environment to install mamba and snakemake into

```
conda env create --name mamba -c conda-forge python3=3.8 mamba
```

2. Activate the new conda environment and install snakemake using mamba

```
conda activate mamba
mamba install -c bioconda -c conda-forge snakemake=7.3.8
```

3. Confirm that snakeake has installed.

```
snakemake --version
```

## Setup

1. Move data into the data/ subdirectory in a folder named, e.g., Batch_01

2. Go into the data/ subdirectory and use the `prep_sample_sheet.sh` script to prepare a sample sheet

```
cd data/
../workflow/scripts/prep_sample_sheet.sh 
```

3. Move the output file into the config folder and confirm in config.yaml file that the samples option is correctly set to that file name

```
mv samples_minimal.tsv ../congfig/samples_minimal.tsv
```
## Execute

1. (Optional) Check that snakemake is correctly interpretting your sample spreadsheet by executing a dryrun or one of the commands in the notes above

```
snakemake --dry-run
```

2. In the base snakemake/project directory, use the `run_snakemake.sh` file to submit the snakemake scheduler as a job to slurm 

```
sbatch run_snakeake.sh
```

alternatively, you can run an interactive and execute the contents of the `run_snakemake.sh` file if you are running into challenges with slurm

3. wait for data :)
(hopefully)

