# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.
import itertools
import os 
import glob
import sys
import pandas as pd

configfile: "config/config.yaml"
include: "workflow/rules/common.smk"
#include: "workflow/rules/metagenomics.smk"
# include: "workflow/rules/bin_assembly.smk"
#include: "workflow/rules/single_sample_binning.smk"
include: "workflow/rules/trim.smk"
include: "workflow/rules/read_qc.smk"
include: "workflow/rules/deconvolute.smk"
include: "workflow/rules/metaphlan.smk"
include: "workflow/rules/assemble.smk"
include: "workflow/rules/assembly_qc.smk"
include: "workflow/rules/spades.smk"
include: "workflow/rules/bin_metabat2.smk"

rule all:
    input:
        get_rules,
        "results/bowtie_out/B16_LyPMA/B16_LyPMA.fastp_bowtie.r1.fastq",
        "results/bowtie_out/B16_LyPMA/B16_LyPMA.fastp_bowtie.r2.fastq"



# Make report for snakemake. 

report: "workflow/report/workflow.rst"
