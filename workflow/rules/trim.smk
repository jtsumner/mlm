import glob
import pandas as pd
from snakemake.utils import validate


rule fastp_pe:
    input:
        r1 = get_r1,
        r2 = get_r2
    output:
        r1Filtered = "results/fastp_out/{sample}/{sample}.fastp.r1.fastq.gz",
        r2Filtered = "results/fastp_out/{sample}/{sample}.fastp.r2.fastq.gz",
        json = "results/fastp_out/{sample}/{sample}_fastp.json",
        html = "results/fastp_out/{sample}/{sample}_fastp.html"
    conda:
        "../envs/seq_processing.yml"
    threads: 16
    shell: 
        "fastp -i {input.r1} -I {input.r2} --out1 {output.r1Filtered} --out2 {output.r2Filtered} --detect_adapter_for_pe --thread {threads} --length_required 50 -j {output.json} -h {output.html} -V"


rule fastqc:
    input: 
        "results/fastp_out/{sample}/{sample}.fastp.r1.fastq.gz",
        "results/fastp_out/{sample}/{sample}.fastp.r2.fastq.gz"
    output:
        "results/fastqc_out/fastp_qc/{sample}.fastp.r1_fastqc.html",
        "results/fastqc_out/fastp_qc/{sample}.fastp.r2_fastqc.html"
    params:
        outDir = "results/fastqc_out/fastp_qc/"
    threads: 12
    shell:
        "module load fastqc/0.11.5 ; fastqc -t {threads} {input} --outdir {params.outDir}"

