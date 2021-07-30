from snakemake.utils import validate
import pandas as pd

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity

samples = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)

def get_read_path(wildcards):
    tmp = samples.loc[wildcards.sample, ["sample", "dataset"]].dropna()
    path = "../data/{}/{}".format(tmp["dataset"], tmp["sample"])
    return path

def get_r1(wildcards):
    tmp1 = get_read_path(wildcards)
    return "{}_R1_001.fastq.gz".format(tmp1)
 
def get_r2(wildcards):
    tmp2 = get_read_path(wildcards)
    return "{}_R2_001.fastq.gz".format(tmp2)

