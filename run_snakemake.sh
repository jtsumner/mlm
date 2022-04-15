#! /bin/bash
#SBATCH -A p31588
#SBATCH --job-name="scheduler"
#SBATCH -t 04:00:00
#SBATCH -N 1
#SBATCH -n 5
#SBATCH -p short
#SBATCH --mem=20gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jacksumner2026@u.northwestern.edu
# module load anaconda3
# source activate snakemake

module purge all
module load anaconda3
source activate mamba

cd $SLURM_SUBMIT_DIR


echo Starting snakemake on the cluster

snakemake --profile simple

