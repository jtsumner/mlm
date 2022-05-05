import glob
import pandas as pd
from snakemake.utils import validate


rule fastp_pe:
    input:
        r1 = get_r1,
        r2 = get_r2
    output:
        r1Filtered = "results/{dataset}/filtered/{sample}.filtered.R1.fastq.gz",
        r2Filtered = "results/{dataset}/filtered/{sample}.filtered.R2.fastq.gz",
        json = "results/{dataset}/filtered/{sample}_fastp.json",
        html = "results/{dataset}/filtered/{sample}_fastp.html"
    conda:
        "../envs/seq_processing.yml"
    threads: 16
    shell: 
        "fastp -i {input.r1} -I {input.r2} --out1 {output.r1Filtered} --out2 {output.r2Filtered} --detect_adapter_for_pe --thread {threads} --length_required 50 -j {output.json} -h {output.html} -V"


rule fastqc:
    input: 
        "results/{dataset}/filtered/{sample}.filtered.R1.fastq.gz",
        "results/{dataset}/filtered/{sample}.filtered.R2.fastq.gz"
    output:
        "results/{dataset}/filtered/fastqc/{sample}.filtered.R1_fastqc.html",
        "results/{dataset}/filtered/fastqc/{sample}.filtered.R2_fastqc.html"
    params:
        outDir = "results/{dataset}/filtered/fastqc"
    threads: 12
    shell:
        "module load fastqc/0.11.5 ; fastqc -t {threads} {input} --outdir {params.outDir}"

