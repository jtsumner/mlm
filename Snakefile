# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.
import itertools
import os 
import glob
import sys
import pandas as pd

configfile: "config/config.yaml"
include: "workflow/rules/common.smk"
include: "workflow/rules/metagenomics.smk"
# include: "workflow/rules/bin_assembly.smk"
include: "workflow/rules/single_sample_binning.smk"

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


rule all:
    input:
        get_rules


# Make report for snakemake. 

report: "report/workflow.rst"
