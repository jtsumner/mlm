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

rule all:
    input:
        get_rules,
        "results/fastqc_out/fastp_qc/DNA_B01_18/DNA_B01_18.fastp.r1_fastqc.html",
        "results/fastqc_out/bbduk_qc/DNA_B01_18/DNA_B01_18.bbduk.r1_fastqc.html",
        "results/fastqc_out/bowtie_qc/DNA_B01_18/DNA_B01_18.bowtie.r1_fastqc.html",
        "results/fastqc_out/bbmerge_qc/DNA_B01_18/DNA_B01_18.bbmerge_fastqc.html",
        "results/fastqc_out/raw_qc/DNA_B01_18/DNA_B01_18.raw.r1_fastqc.html",
        "results/fastqc_out/multiqc_report.html"

# Make report for snakemake. 
report: "workflow/report/workflow.rst"

