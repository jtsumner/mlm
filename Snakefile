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

rule all:
    input:
        # expand("../results/filtered/{dataset}/{sample}.filtered.R1.fastq.gz",
        #     sample=samples["sample"], dataset=samples["dataset"]),
        # expand("../results/filtered/{dataset}/{sample}.filtered.R2.fastq.gz", 
        #     sample=samples["sample"], dataset=samples["dataset"])
        "results/allDatasets/metaphlan/abundance_heatmap_species.allDatasets.png",
        "results/allDatasets/metaphlan/abundance_heatmap_genus.allDatasets.png" # hclust genus 
        # #"results/allDatasets/coassembly/quast/report.html",
        # expand("results/{dataset}/assembly/quast/{sample}_quast/report.html",
        #     zip, sample=samples["sample"], dataset=samples["dataset"]),
        # "results/allDatasets/metaphlan/unifrac_matrix.allDatasets.txt",
        # "results/allDatasets/single_sample_assemblies/multiqc_stats/report.html"
        # "results/allDatasets/single_sample_assemblies/allSamples.megahit_g1000.fa" # merged clean monoassemblies
        # "results/allDatasets/single_sample_assemblies/quast/monoassemblies_quast/report.html"
        # # "results/allDatasets/single_sample_assemblies/metabat2/bins/bin.1.fa", # binned_contigs 
        # "results/allDatasets/single_sample_assemblies/metabat2/plots/bin_qa_plot.png"
        # # expand("results/{dataset}/assembly/megahit_g1000/metabat2/checkm/{sample}/{sample}_checkm_output.txt",
        # #     zip, sample=samples["sample"], dataset=samples["dataset"])


# Make report for snakemake. 

#report: "report/workflow.rst"
