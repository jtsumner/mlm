#! /bin/bash
#SBATCH -A p31288
#SBATCH --job-name="report"
#SBATCH -t 01:00:00
#SBATCH -n 2
#SBATCH -p short
#SBATCH --mem-per-cpu=3gb

module load anaconda3
source activate snakemake

# Must be in microbiome-snakemake/workflow/ directory to execute
cd $SLURM_SUBMIT_DIR

#--max-jobs-per-second 5 --max-status-checks-per-second 5 
snakemake --report report.html
