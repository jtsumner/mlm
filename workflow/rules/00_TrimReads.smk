import glob
import pandas as pd
from snakemake.utils import validate

############################
### PART 1A: FASTP TRIM  ###
############################

rule fastp_pe:
    input:
        r1 = get_r1,
        r2 = get_r2
    output:
        r1_filtered = "results/fastp_out/{sample}/{sample}.fastp.r1.fastq.gz",
        r2_filtered = "results/fastp_out/{sample}/{sample}.fastp.r2.fastq.gz",
        json = "results/fastp_out/{sample}/{sample}_fastp.json",
        html = "results/fastp_out/{sample}/{sample}_fastp.html"
    conda:
        "../envs/seq_processing.yml"
    threads: 12
    resources:
        mem="20G"
    shell: 
        """
        fastp \
            -i {input.r1} \
            -I {input.r2} \
            --out1 {output.r1_filtered} \
            --out2 {output.r2_filtered} \
            --detect_adapter_for_pe \
            --trim_poly_g \
            --trim_poly_x \
            --dedup \
            --thread {threads} \
            --length_required 50 \
            -j {output.json} \
            -h {output.html} \
            -V 
        """

############################
###  PART 2: COMPLEXITY  ###
############################

rule complexity_filter:
    """
    For entropy-based filtering of low complexity reads
    Good for host-associated metagenomics
    """
    input:
        r1 = get_trimmed_read1,
        r2 = get_trimmed_read2
    output:
        r1 = "results/bbduk_out/{sample}/{sample}.bbduk.r1.fastq.gz",
        r2 = "results/bbduk_out/{sample}/{sample}.bbduk.r2.fastq.gz"
    params:
        entropy = 0.5
    threads: 5
    resources:
        mem="10G"
    shell: 
        """
        module load BBMap
        bbduk.sh \
            in={input.r1} \
            in2={input.r2} \
            out={output.r1} \
            out2={output.r2} \
            entropy={params.entropy} \
            entropywindow=50 \
            entropyk=5 \
            threads={threads}
        """

############################
## PART 3: DECONTAMINATE  ##
############################

rule qc_filter:
    """
    Performs host read filtering on paired end data using Bowtie and Samtools
    """
    input:
        r1 = get_penultimate_read1,
        r2 = get_penultimate_read2
    output:
        r1_clean = "results/bowtie_out/{sample}/{sample}.bowtie.r1.fastq.gz",
        r2_clean = "results/bowtie_out/{sample}/{sample}.bowtie.r2.fastq.gz",
        sorted_bam = "results/bowtie_out/{sample}/{sample}.mapped.sorted.bam",
        flagstat = "results/bowtie_out/{sample}/{sample}.flagstat.tsv"
    params:
        filter_db = "resources/bowtie_human/GRCh38_noalt_as/GRCh38_noalt_as",
        #filter_db = "/projects/b1042/HartmannLab/jack/mlm/resources/bowtie_mouse/GRCm39/GRCm39",
    threads: 15
    resources:
        mem="25G",
        time="02:00:00"
    shell:
        """
        module load bowtie2
        module load samtools/1.10.1

        bowtie2 -p {threads} -x {params.filter_db} --very-sensitive -1 {input.r1} -2 {input.r2}| \
        samtools view -bS -@ {threads}| \
        samtools sort -@ {threads} -n -o {output.sorted_bam}

        samtools fastq -1 {output.r1_clean} -2 {output.r2_clean} -@ {threads} -f 12 -F 256 {output.sorted_bam}

        samtools flagstat -@ {threads} -O tsv {output.sorted_bam} > {output.flagstat}
        """

############################
### PART 4: MERGE READS  ###
############################

rule merge_reads:
    """
    For merging reads PE reads to improve alignment rates 
    """
    input:
        r1 = get_final_read1,
        r2 = get_final_read2
    output:
        r3 = "results/bbmerge_out/{sample}/{sample}.bbmerge.fastq.gz"
    threads: 5
    resources:
        mem="10G"
    shell:
        """
        module load BBMap
        bbmerge.sh \
            in1={input.r1} \
            in2={input.r2} \
            out={output.r3} \
            strict=t \
            k=60 \
            mininsert=70 \
            threads={threads}
        """