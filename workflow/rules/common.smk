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
#     path = "./data  /{}/{}".format(dt, samp)
#     return path

# def get_r1(wildcards):
#     tmp1 = get_read_path(wildcards)
#     return "{}_R1_001.fastq.gz".format(tmp1)
 
# def get_r2(wildcards):
#     tmp2 = get_read_path(wildcards)
#     return "{}_R2_001.fastq.gz".format(tmp2)

def get_rules(wildcards):
    all_rules = []
    if config["METAPHLAN"]:
        all_rules.append("results/allDatasets/metaphlan/abundance_heatmap_species.allDatasets.png")
        all_rules.append("results/allDatasets/metaphlan/abundance_heatmap_genus.allDatasets.png")
        all_rules.append("results/allDatasets/metaphlan/unifrac_matrix.allDatasets.txt")
    if config["ASSEMBLE"]:
        megahit_results = expand("results/{dataset}/assembly/quast/{sample}_quast/report.html", zip, sample=samples["sample"], dataset=samples["dataset"])
        all_rules = all_rules + megahit_results
        all_rules.append("results/allDatasets/single_sample_assemblies/multiqc_stats/report.html")
    if config["METABAT2"]:
        metabat2_results = expand("results/{dataset}/assembly/megahit_g1000/metabat2/{sample}/bins/bin.1.fa", zip, sample=samples["sample"], dataset=samples["dataset"])
        all_rules = all_rules + metabat2_results
    return all_rules


def get_read_path(wildcards):
    return samples.loc[wildcards.sample, ["sample", "dataset"]].dropna()

def get_r1(wildcards):
    tmp = get_read_path(wildcards)
    if tmp["dataset"] == "Batch_04":
        return "data/" + tmp["dataset"] + "/" + tmp["sample"] + "_R1.fastq.gz"
    return "data/" + tmp["dataset"] + "/" + tmp["sample"] + "_R1_001.fastq.gz"

def get_r2(wildcards):
    tmp = get_read_path(wildcards)
    if tmp["dataset"] == "Batch_04":
        return "data/" + tmp["dataset"] + "/" + tmp["sample"] + "_R2.fastq.gz"
    return "data/" + tmp["dataset"] + "/" + tmp["sample"] + "_R2_001.fastq.gz"
