# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.
import itertools
import os 
import glob
import sys
import pandas as pd

configfile: "config/config.yaml"
include: "workflow/rules/common.smk"
include: "workflow/rules/00_TrimReads.smk"
include: "workflow/rules/01_DecontaminateReads.smk"
include: "workflow/rules/02_TaxonomicAnalysis.smk"
include: "workflow/rules/04_ViralAnalysis.smk"
include: "workflow/rules/05_AssemblyAnalysis.smk"
include: "workflow/rules/06_SingleGenomeBins.smk"


rule all:
    input:
        get_rules,
        #"results/humann_out/DNA_B03_02/DNA_B03_02_genefamilies.tsv"
        #expand("results/humann_out/{sample}/{sample}_genefamilies.tsv", sample=samples["sample"]), #added
        # expand("results/metaphlan_bbmerge_out/{sample}/{sample}.metaphlan_profile.txt",sample=samples["sample"])# added
        #"results/metaphlan_bowtie_out/merged_metaphlan_profile.tsv",

# Make report for snakemake. 
report: "workflow/report/workflow.rst"

