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
        #get_rules,
        "results/metaphlan_bowtie_out/merged_metaphlan_profile.tsv",
        "results/spades_parsed/B16_LyPMA/B16_LyPMA.fa"
        #expand("results/prokka_out/{sample}/{sample}.tsv", sample=samples["sample"]), # Prokka
        #expand("results/checkm_out/{sample}/{sample}_checkm_output.txt", sample=samples["sample"]), # checkm
        #"results/negative_db/negative_controls.contigs.fa", # make negative control datbase
        #expand("results/negative_out/{sample}/{sample}.clean.r1.fastq.gz", sample=samples["sample"])
        #"results/vcontact2_data/vcontact2_output/genome_by_genome_overview.csv"
# Make report for snakemake. 
report: "workflow/report/workflow.rst"

