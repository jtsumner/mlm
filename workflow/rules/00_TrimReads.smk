import glob
import pandas as pd
from snakemake.utils import validate
import os.path
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
            --thread {threads} \
            --length_required 50 \
            -j {output.json} \
            -h {output.html} \
            -V 
        """

rule fastp4amp:
    """
    This is a special trimming rule used to quality trim and merge 
    reads from 16S amplicon sequencing. The remainder of my 16S pipeline 
    is not integrated into MLM (see MGX2AMP repo).

    After this step, prepare manifest and sample metadata files,
    acquire the SILVA DB for 16S sequences, and run QIIME2

    To use: 
        1. Comment out get_rules in rule all.
        2. Replace with expand(
            "results/AMP_trimmed/{sample}_fastp-merged.fq.gz", 
            sample=samples["sample"])
        3. Prep sample sheet and execute
    """
    input:
        r1 = get_r1,
        r2 = get_r2
    output:
        r1_filtered = "results/AMP_trimmed/{sample}/{sample}.fastp.r1.fastq.gz",
        r2_filtered = "results/AMP_trimmed/{sample}/{sample}.fastp.r2.fastq.gz",
        json = "results/AMP_trimmed/{sample}/{sample}_fastp.json",
        html = "results/AMP_trimmed/{sample}/{sample}_fastp.html"
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
            --thread {threads} \
            --length_required 100 \
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
        entropy = 0.3
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
            entropywindow=100 \
            entropyk=10 \
            threads={threads}
        """

############################
## PART 3: DECONTAMINATE  ##
############################

rule host_decontamination:
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

        bowtie2 -p {threads} -x {params.filter_db} --local -1 {input.r1} -2 {input.r2}| \
        samtools view -bS -@ {threads}| \
        samtools sort -@ {threads} -n -o {output.sorted_bam}

        samtools fastq -1 {output.r1_clean} -2 {output.r2_clean} -@ {threads} -f 12 -F 256 {output.sorted_bam}

        samtools flagstat -@ {threads} -O tsv {output.sorted_bam} > {output.flagstat}
        """

rule flagstat_summarize:
    input:
        expand("results/bowtie_out/{sample}/{sample}.flagstat.tsv",
        sample=samples["sample"])
    output:
        "results/bowtie_out/flagstat_summary.txt"
    shell:
        """
        cd results/bowtie_out/
        for i in $(ls -d *) ; do sed -e "s/^/$i\t/" $i/*.flagstat.tsv >> flagstat_summary.txt ; done
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
            mininsert=90 \
            threads={threads}
        """

###############################
###  PART 5A: QC NONPAREIL  ###
###############################

rule nonpareil:
    """
    Takes fastq.gz files from bowtie2 output 
    Decompresses reads and pipes stdout to a tmp file in nonpareil folder
    Kmer based nonpareil executed on tmp file and them tmp file is deleted
    """
    input:
        r1_clean = get_final_read1,
        r2_clean = get_final_read2,
    output:
        tmp_fastq = temp("results/nonpareil_out/{sample}/{sample}.tmp"),
        summary = "results/nonpareil_out/{sample}/{sample}.npo",
        out_dir = directory("results/nonpareil_out/{sample}/")
    threads: 16
    resources:
        mem="25G"
    shell:
        """
        module load nonpareil/3.4.1
        gunzip -c {input.r1_clean} > {output.tmp_fastq}
        nonpareil -s {output.tmp_fastq} -T kmer -f fastq -b {output.out_dir}/{wildcards.sample} -t {threads}
        """
        
use rule nonpareil as nonpareil_merged with:
    input:
        r1_clean = get_final_merged_read,
        r2_clean = get_final_read2
    output:
        tmp_fastq = temp("results/nonpareil.merged_out/{sample}/{sample}.tmp"),
        summary = "results/nonpareil.merged_out/{sample}/{sample}.npo",
        out_dir = directory("results/nonpareil.merged_out/{sample}/")

###############################
###    PART 5B: QC FASTQC   ###
###############################

rule fastqc_fastp:
    input: 
        get_trimmed_read1,
        get_trimmed_read2
    output:
        "results/fastqc_out/fastp_qc/{sample}/{sample}.fastp.r1_fastqc.html",
        "results/fastqc_out/fastp_qc/{sample}/{sample}.fastp.r2_fastqc.html"
    params:
        out_dir = "results/fastqc_out/fastp_qc/{sample}"
    threads: 4
    resources:
        time = "00:30:00"
    shell:
        """
        module load fastqc/0.11.5
        fastqc -t {threads} {input} --outdir {params.out_dir}
        """ 

use rule fastqc_fastp as fastqc_bbduk with:
    input: 
        get_complex_read1,
        get_complex_read2
    output:
        "results/fastqc_out/bbduk_qc/{sample}/{sample}.bbduk.r1_fastqc.html",
        "results/fastqc_out/bbduk_qc/{sample}/{sample}.bbduk.r2_fastqc.html"
    params:
        out_dir = "results/fastqc_out/bbduk_qc/{sample}"

use rule fastqc_fastp as fastqc_bowtie with:
    input: 
        get_decontaminated_read1,
        get_decontaminated_read2
    output:
        "results/fastqc_out/bowtie_qc/{sample}/{sample}.bowtie.r1_fastqc.html",
        "results/fastqc_out/bowtie_qc/{sample}/{sample}.bowtie.r2_fastqc.html"
    params:
        out_dir = "results/fastqc_out/bowtie_qc/{sample}"

use rule fastqc_fastp as fastqc_bbmerge with:
    input: 
        get_merged_reads
    output:
        "results/fastqc_out/bbmerge_qc/{sample}/{sample}.bbmerge_fastqc.html"
    params:
        out_dir = "results/fastqc_out/bbmerge_qc/{sample}"

rule fastqc_raw:
    input: 
        r1 = get_r1,
        r2 = get_r2
    output:
        html1 = "results/fastqc_out/raw_qc/{sample}/{sample}.raw.r1_fastqc.html",
        html2 = "results/fastqc_out/raw_qc/{sample}/{sample}.raw.r2_fastqc.html",
        zip1 = "results/fastqc_out/raw_qc/{sample}/{sample}.raw.r1_fastqc.zip",
        zip2 = "results/fastqc_out/raw_qc/{sample}/{sample}.raw.r2_fastqc.zip",
        outdir = directory("results/fastqc_out/raw_qc/{sample}/")
    params:
        html1 = lambda wildcards, output, input: output.outdir + "/" + os.path.basename(input.r1).split(".fastq.gz")[0] + "_fastqc.html",
        html2 = lambda wildcards, output, input: output.outdir + "/" + os.path.basename(input.r2).split(".fastq.gz")[0] + "_fastqc.html",
        zip1 = lambda wildcards, output, input: output.outdir + "/" + os.path.basename(input.r1).split(".fastq.gz")[0] + "_fastqc.zip",
        zip2 = lambda wildcards, output, input: output.outdir + "/" + os.path.basename(input.r2).split(".fastq.gz")[0] + "_fastqc.zip",
    threads: 4
    resources:
        time = "00:30:00"
    shell:
        """
        module load fastqc/0.11.5
        mkdir -p {output.outdir}
        fastqc -t {threads} {input.r1} {input.r2} --outdir {output.outdir}
        mv {params.html1} {output.html1}
        mv {params.html2} {output.html2}
        mv {params.zip1} {output.zip1}
        mv {params.zip2} {output.zip2}
        """

rule fastqc_multiqc:
    input:
        get_multiqc_input()
    output:
        multiqc_report = "results/fastqc_out/multiqc_report.html"
    params:
        out_dir="results/fastqc_out"
    shell:
        """
        module load multiqc
        multiqc --outdir {params.out_dir} --dirs --dirs-depth 2 results/fastqc_out/ -f
        """

###############################
##   PART 5B: CONTAMINANTS   ##
###############################

def get_negative_reads1(wildcards):
    r1 = get_trimmed_read1(wildcards)
    return expand(r1, sample=samples[samples["sample"].isin(negative_controls)]["sample"])

def get_negative_reads2(wildcards):
    r2 = get_trimmed_read2(wildcards)
    return expand(r2, sample=samples[samples["sample"].isin(negative_controls)]["sample"])

rule get_control_assemblies:
    """
    TAKES trimmed reads from negative_controls (see config)
    RETURNS coassembled and bowtie indexed database from contaminant removal

    Functionality in params.r1 and params.r2 to take lists from input and convert to 
    comma-seperated strings for direct input to megahit, bypassing the need to 
    concatenate files manually. Manually confirmed r1 and r2 are in same order
    as expected
    """
    input:
        r1 = get_negative_reads1,
        r2 = get_negative_reads2
    output:
        out_dir = directory("results/negative_db"),
        negative_scaffolds = "results/negative_db/negative_controls.contigs.fa"
    params:
        r1 = ",".join(get_negative_reads1(samples)),
        r2 = ",".join(get_negative_reads2(samples))
    resources:
        mem="130g",
        time="02:00:00"
    threads: 40
    shell:
        """
        module load megahit/1.0.6.1 bowtie2/2.4.5 

        megahit -t {threads} \
            -m 120e9 \
            -1 {params.r1} \
            -2 {params.r2} \
            --out-prefix negative_controls -o {output.out_dir}/assemble \
            --presets meta \
            --min-contig-len 1000
        mv {output.out_dir}/assemble/negative_controls.contigs.fa {output.negative_scaffolds} 
        bowtie2-build {output.negative_scaffolds} {output.out_dir}/negative_db --thread {threads} --seed 12345
        """

#Performs host read filtering on paired end data using Bowtie and Samtools 
use rule host_decontamination as negative_decontamination with:
    input:
        r1 = get_penultimate_read1,
        r2 = get_penultimate_read2,
        negative_scaffolds = "results/negative_db/negative_controls.contigs.fa"
    output:
        r1_clean = "results/negative_out/{sample}/{sample}.clean.r1.fastq.gz",
        r2_clean = "results/negative_out/{sample}/{sample}.clean.r2.fastq.gz",
        sorted_bam = "results/negative_out/{sample}/{sample}.mapped.sorted.bam",
        flagstat = "results/negative_out/{sample}/{sample}.flagstat.tsv"
    params:
        filter_db = "results/negative_db/negative_db"
    threads: 15
    resources:
        mem="25G",
        time="02:00:00"
