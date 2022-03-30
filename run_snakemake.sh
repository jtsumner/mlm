#! /bin/bash
#SBATCH -A p31588
#SBATCH --job-name="scheduler"
#SBATCH -t 04:00:00
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -p short
#SBATCH --mem=5gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jacksumner2026@u.northwestern.edu
# module load anaconda3
# source activate snakemake

module purge all
module load anaconda3
#source activate mamba
#source activate updated_snakemake
#mamba activate updated_snakemake
source activate snakemake

# Must be in microbiome-snakemake/workflow/ directory to execute
cd $SLURM_SUBMIT_DIR

#--max-jobs-per-second 5 --max-status-checks-per-second 5 
mkdir -p logs_slurm
# snakemake --verbose --use-conda --cluster-config cluster.yaml --max-jobs-per-second 1 --max-status-checks-per-second 1 -j 100 --cluster "sbatch -A {cluster.allocation} -p {cluster.partition} -t {cluster.time} --mem={cluster.mem} -N {cluster.nodes} -n {cluster.cpus} -o {cluster.output} -e {cluster.error} --mail-type={cluster.email_type} --mail-user={cluster.email} --job-name={cluster.jobname}"


echo Starting snakemake on the cluster
snakemake --verbose \
    --use-conda \
    --cluster-config config/cluster.yaml \
    --max-jobs-per-second 1 \
    --max-status-checks-per-second 1 \
    -j 100 \
    --cluster "sbatch -A {cluster.allocation} -p {cluster.partition} -t {cluster.time} --mem={cluster.mem} -N {cluster.nodes} -n {cluster.cpus} -o {cluster.output} -e {cluster.error} --mail-type={cluster.email_type} --mail-user={cluster.email} --job-name={cluster.jobname}"
