import glob
import pandas as pd
from snakemake.utils import validate
import os.path


############################
### PART 1A: FASTP TRIM  ###
############################

rule rna_fastp:
    input:
        r1 = get_r1,
        r2 = get_r2
    output:
        r1_filtered = "results/rna/fastp_out/{sample}/{sample}.fastp.r1.fastq.gz",
        json = "results/rna/fastp_out/{sample}/{sample}_fastp.json",
        html = "results/rna/fastp_out/{sample}/{sample}_fastp.html"
    conda:
        "../envs/seq_processing.yml"
    threads: 12
    resources:
        mem="20G"
    shell: 
        """
        fastp \
            -i {input.r1} \
            --out {output.r1_filtered} \
            --detect_adapter_for_pe \
            --thread {threads} \
            -j {output.json} \
            -h {output.html} \
            -V 
        """


rule rna_host_decontamination:
    """
    Performs host read filtering on paired end data using Bowtie and Samtools
    """
    input:
        r1 = "results/rna/fastp_out/{sample}/{sample}.fastp.r1.fastq.gz"
    output:
        r1_clean = "results/rna/bowtie_out/{sample}/{sample}.bowtie.r1.fastq.gz",
        sorted_bam = "results/rna/bowtie_out/{sample}/{sample}.mapped.sorted.bam",
        flagstat = "results/rna/bowtie_out/{sample}/{sample}.flagstat.tsv"
    params:
        filter_db = "resources/bowtie_human/chm13.draft_v1.0_plusY/chm13.draft_v1.0_plusY"
        #filter_db = "resources/bowtie_human/GRCh38_noalt_as/GRCh38_noalt_as",
        #filter_db = "/projects/b1042/HartmannLab/jack/mlm/resources/bowtie_mouse/GRCm39/GRCm39",
    threads: 15
    resources:
        mem="25G",
        time="02:00:00"
    shell:
        """
        module load bowtie2/2.4.5
        module load samtools/1.10.1

        bowtie2 -p {threads} -x {params.filter_db} --very-sensitive -1 {input.r1}| \
        samtools view -bS -@ {threads}| \
        samtools sort -@ {threads} -n -o {output.sorted_bam}

        samtools fastq -1 {output.r1_clean} -@ {threads} -f 12 -F 256 {output.sorted_bam}

        samtools flagstat -@ {threads} -O tsv {output.sorted_bam} > {output.flagstat}
        """


rule rna_kraken2:
    """
    Performs taxnomic classification with Kraken2 

    Outputs a kraken2-style report and metaphlan-style report with a
    script from KrakenTools
    """
    input: 
        r1_clean = "results/rna/bowtie_out/{sample}/{sample}.bowtie.r1.fastq.gz"
    output:
        kraken_out = "results/rna/kraken/{sample}/{sample}_kraken2out.txt",
        kraken_report = "results/rna/kraken/{sample}/{sample}_kraken2report.txt"
    threads: 20
    resources:
        mem="60G",
        time="01:00:00",
        partition="genomics-himem"
    params:
        kraken_db = "/projects/b1188/bmo/mlm/resources/kraken_db" #/software/kraken/database/kraken_db
    shell:
        """
        module load kraken/2
        kraken2 --threads {threads} \
            --output {output.kraken_out} \
            --report {output.kraken_report} \
            --confidence 0.7 \
            --gzip-compressed \
            --minimum-hit-groups 3 \
            --db {params.kraken_db} {input.r1_clean}
        """ 

rule concat_reads:
    """
    For merging reads PE reads to improve alignment rates 
    """
    input:
        r1 = get_final_read1,
        r2 = get_final_read2
    output:
        r3 = "results/concat_reads/{sample}/{sample}.fastq.gz"
    threads: 5
    resources:
        mem="20G"
    shell:
        """
        cat {input.r1} {input.r2} > {output.r3}
        """