import glob
import pandas as pd
from snakemake.utils import validate



rule fastp_pe:
    input:
        r1 = get_r1,
        r2 = get_r2
    output:
        r1Filtered = "../results/filtered/{dataset}/{sample}.filtered.R1.fastq.gz",
        r2Filtered = "../results/filtered/{dataset}/{sample}.filtered.R2.fastq.gz",
        json = "../results/filtered/{dataset}/{sample}_fastp.json",
        html = "../results/filtered/{dataset}/{sample}_fastp.html"
    conda:
        "../envs/old_envs/seq_processing.yml"
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
    threads: 12
    shell:
        "module load fastqc/0.11.5 ; fastqc -t {threads} {input} --outdir {params.outDir}"

rule bwa_map:
    input:
        r1Filtered = "../results/filtered/{dataset}/{sample}.filtered.R1.fastq.gz",
        r2Filtered = "../results/filtered/{dataset}/{sample}.filtered.R2.fastq.gz"
    output:
        sam = temp("../results/bwa/{dataset}/{sample}.mapped.sam"),
        bam = "../results/bwa/{dataset}/{sample}.mapped.bam",
        sortedBam = "../results/bwa/{dataset}/{sample}.mapped.sorted.bam",
        unmappedBam = "../results/bwa/{dataset}/{sample}.unmapped.bam",
        cleanFastQ1 = "../results/bwa/{dataset}/{sample}.clean.R1.fastq",
        cleanFastQ2 = "../results/bwa/{dataset}/{sample}.clean.R2.fastq"
    params:
        genome = "/projects/b1042/HartmannLab/jack/SCRIPT/expPipeline_v1/data/genome/hg38.fa.gz"
    threads: 20
    shell:
        """
        module purge all
        module load bwa/0.7.17
        module load samtools/1.10.1
        module load bedtools/2.29.2
        bwa mem -t {threads} {params.genome} {input.r1Filtered} {input.r2Filtered} > {output.sam}
        samtools view -Subh -o {output.bam} {output.sam}
        samtools sort -o {output.sortedBam} {output.bam}

        samtools view -b -f 12 -F 256 -o {output.unmappedBam} {output.sortedBam}
        bedtools bamtofastq -i {output.unmappedBam} -fq {output.cleanFastQ1} -fq2 {output.cleanFastQ2}
        """

rule metaphlan:
    input:
        cleanFastQ1 = "../results/bwa/{dataset}/{sample}.clean.R1.fastq",
        cleanFastQ2 = "../results/bwa/{dataset}/{sample}.clean.R2.fastq"
    output:

"""
        # samtools view -Subh {output.sam} | samtools sort –o {output.bam} -

rule multiqc:
    input:
        expand("../results/filtered/{dataset}/{sample}.filtered.R1.fastq.gz"),
        expand("../results/filtered/{dataset}/{sample}.filtered.R2.fastq.gz")
    output:
        expand("../results/filtered/{dataset}/fastqc/multiqc_report.html")
    params:
        inDir="../results/filtered/{datasets}/fastqc"

    shell:
        "module load multiqc; multiqc {params.inDir}"
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

"""