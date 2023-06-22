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
#include: "workflow/rules/read_qc.smk"
include: "workflow/rules/01_DecontaminateReads.smk"
include: "workflow/rules/02_TaxonomicAnalysis.smk"
include: "workflow/rules/assemble.smk"
include: "workflow/rules/assembly_qc.smk"
include: "workflow/rules/bin_metabat2.smk"
include: "workflow/rules/04_ViralAnalysis.smk"

rule all:
    input:
        get_rules,
        "results/vcontact2_data/vcontact2_output/genome_by_genome_overview.csv"        #"results/vibrant_output/VIBRANT_parsed_scaffolds/VIBRANT_phages_parsed_scaffolds/parsed_scaffolds.phages_combined.faa",
        #"results/vibrant_output/VIBRANT_parsed_scaffolds/VIBRANT_phages_parsed_scaffolds/parsed_scaffolds.phages_combined.simple.faa"

        #expand("results/vibrant_output/{sample}/VIBRANT_{sample}/VIBRANT_phages_{sample}/{sample}.phages_combined.faa", sample=samples["sample"]),
        #"results/spades_parsed/DNA_B03_27/DNA_B03_27.fa", 
        #"results/vibrant_output/DNA_B03_27/VIBRANT_DNA_B03_27/VIBRANT_phages_DNA_B03_27/DNA_B03_27.phages_combined.faa"
        #expand("results/AMP_trimmed/{sample}_fastp-merged.fq.gz", sample=samples["sample"])        #g

# Make report for snakemake. 
report: "workflow/report/workflow.rst"

