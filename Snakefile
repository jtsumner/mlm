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

rule all:
    input:
        get_rules


# Make report for snakemake. 

report: "report/workflow.rst"
