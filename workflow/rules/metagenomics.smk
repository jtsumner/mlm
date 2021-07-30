import glob

import pandas as pd
from snakemake.utils import validate


validate(config, schema="../schemas/config.schema.yaml")

samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"sample_name": str})
    .set_index("sample_name", drop=False)
    .sort_index()
)

def get_fastp_input(wildcards):
    pass


rule fastp_pe:
    input:
    output:
    log:
    conda:
    params:
    threads:
    shell:

rule fastqc
    input:
    output:
    log:
    conda:
    params:
    threads:
    shell:

rule multiqc:
    input:
    output:
    log:
    conda:
    params:
    threads:
    shell:

### deconvolution

rule get_genome:
    input:
    output:
    log:
    conda:
    params:
    threads:
    shell:

rule bwa_index:
    input:
    output:
    log:
    conda:
    params:
    threads:
    shell:

rule bwa_unmapped:
    input:
    output:
    log:
    conda:
    params:
    threads:
    shell:

rule get_kneaddata_db:
    input:
    output:
    log:
    conda:
    params:
    threads:
    shell:

rule kneaddata:
    input:
    output:
    log:
    conda:
    params:
    threads:
    shell:

rule metaphlan3:
    input:
    output:
    log:
    conda:
    params:
    threads:
    shell:

