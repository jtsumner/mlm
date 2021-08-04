from snakemake.utils import validate
import pandas as pd

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity

samples = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
samples.index.names = ["sample"]

# Call and define initial fastp analysis
# def get_read_path(wildcards):
#     tmp = samples.loc[wildcards.sample, ["sample", "dataset"]].dropna()
#     samp = tmp["sample"]
#     dt = tmp["dataset"]
#     path = "../data/{}/{}".format(dt, samp)
#     return path

# def get_r1(wildcards):
#     tmp1 = get_read_path(wildcards)
#     return "{}_R1_001.fastq.gz".format(tmp1)
 
# def get_r2(wildcards):
#     tmp2 = get_read_path(wildcards)
#     return "{}_R2_001.fastq.gz".format(tmp2)

def get_read_path(wildcards):
    return samples.loc[wildcards.sample, ["sample", "dataset"]].dropna()

def get_r1(wildcards):
    tmp = get_read_path(wildcards)
    return "../data/" + tmp["dataset"] + "/" + tmp["sample"] + "_R1_001.fastq.gz"

def get_r2(wildcards):
    tmp = get_read_path(wildcards)
    return "../data/" + tmp["dataset"] + "/" + tmp["sample"] + "_R2_001.fastq.gz"
