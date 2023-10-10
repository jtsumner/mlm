#! /bin/bash
#SBATCH -A b1042
#SBATCH -p genomics-himem
#SBATCH --job-name="scheduler"
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --mem=15G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jts.quest.notifications@gmail.com
#SBATCH --output="slurm_mlm.out"
# module load anaconda3
# source activate snakemake

cd $SLURM_SUBMIT_DIR

# Annotating the output file
START_TIME=$(date)
cat workflow/examples/ascii_art_flowers.txt
echo "
NEW SNAKEMAKE EXECUTION :)
Job Details
Job ID: ${SLURM_JOB_ID}
Start Time: ${START_TIME}

Loading conda...
"

# Load Conda Environment with Snakemake
module purge all
module load mamba
#mamba init
mamba activate snakemake
#source activate snakemamba

# Execute snakemake
echo "Starting snakemake on cluster..."
snakemake --profile simple --prioritize renorm_humann renorm_humann_path metaphlan_bowtie_species_abundance regroup_humann

snakemake --forceall --rulegraph | dot -Tpdf > results/dag.pdf
# Annotating the output file
END_TIME=$(date)
echo "
ENDING SNAKEMAKE EXECUTION
Job Details
Job ID: ${SLURM_JOB_ID}
Start Time: ${END_TIME}

Bye-bye :)
"
