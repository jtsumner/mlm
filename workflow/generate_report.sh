#! /bin/bash
#SBATCH -A p31288
#SBATCH --job-name="report"
#SBATCH -t 01:00:00
#SBATCH -n 2
#SBATCH -p short
#SBATCH --mem-per-cpu=3gb

module purge all
module load python-miniconda3/4.12.0
source activate snakemamba

# Must be in microbiome-snakemake/workflow/ directory to execute
cd $SLURM_SUBMIT_DIR
#sed -i 's/"starttime": null/"starttime": 0/g' *
#--max-jobs-per-second 5 --max-status-checks-per-second 5 
snakemake --report report.html
