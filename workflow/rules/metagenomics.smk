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
        r1 = "../data/{dataset}/{sample}_R1_001.fastq.gz",
        r2 = "../data/{dataset}/{sample}_R2_001.fastq.gz"
    output:
        r1Filtered = "../results/filtered/{dataset}/{sample}.filtered.R1.fastq.gz",
        r2Filtered = "../results/filtered/{dataset}/{sample}.filtered.R2.fastq.gz",
        json = "../results/filtered/{dataset}/{sample}_fastp.json",
        html = "../results/filtered/{dataset}/{sample}_fastp.html"
    threads: 16
    shell: 
        "fastp -i {input.r1} -I {input.r2} --out1 {output.r1Filtered} --out2 {output.r2Filtered} --detect_adapter_for_pe --thread {threads} --length_required 50 -j {output.json} -h {output.html} -V"

rule fastqc:
    input: 
        "../results/filtered/{dataset}/{sample}.filtered.R1.fastq.gz",
        "../results/filtered/{dataset}/{sample}.filtered.R2.fastq.gz"

    output:
        "../results/filtered/{dataset}/fastqc/{sample}.filtered.R1_fastqc.html",
        "../results/filtered/{dataset}/fastqc/{sample}.filtered.R2_fastqc.html"

    params:
        outDir = "../results/filtered/{dataset}/fastqc"
#    wildcard_constraints:
#        reads="[R]1|2"
    threads: 12
    shell:
        "module load fastqc/0.11.5; fastqc -t {threads} {input} --outdir {params.outDir}"

rule multiqc:
    input:
        expand("../results/filtered/{dataset}/fastqc/{sample}.filtered.{read}_fastqc.html", sample=samples["sample_name"], read=READS,             dataset=DATASETS, 
)
    output:
        "../results/filtered/{dataset}/fastqc/multiqc_report.html"
    params:
        inDir="../results/filtered/{datasets}/fastqc"

    shell:
        "module load multiqc; multiqc {params.inDir}"

### deconvolution
"""
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

"""
