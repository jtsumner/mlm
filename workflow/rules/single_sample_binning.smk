
import glob
import pandas as pd
from snakemake.utils import validate

### For single sample assemblies 

rule SS_index_contigs:
    input:
        "../results/{dataset}/assembly/megahit_g1000/{sample}.megahit_g1000.fa"
    output:
        "../results/{dataset}/assembly/megahit_g1000/{sample}.megahit_g1000.bwt"
    params:
        prefix = "../results/{dataset}/assembly/megahit_g1000/{sample}.megahit_g1000"
    shell:
        """
        module load bwa/0.7.17
        module load samtools/1.10.1

        bwa index -p {params.prefix} {input}
        """


rule SS_map_reads_to_contigs:
    input:
        cleanFastQ1 = "../results/{dataset}/bwa/{sample}.clean.R1.fastq",
        cleanFastQ2 = "../results/{dataset}/bwa/{sample}.clean.R2.fastq",
        index = "../results/{dataset}/assembly/megahit_g1000/{sample}.megahit_g1000.bwt"
    output:
        sortedBam = "../results/{dataset}/assembly/megahit_g1000/mapped_reads/{sample}.mapped.sorted.bam"
    params:
        genome = "../results/{dataset}/assembly/megahit_g1000/{sample}.megahit_g1000",
        sam = "../results/{dataset}/assembly/megahit_g1000/mapped_reads/{sample}.mapped.sam",
        bam = "../results/{dataset}/assembly/megahit_g1000/mapped_reads/{sample}.mapped.bam",
    threads: 20
    shell:
        """
        module purge all
        module load bwa/0.7.17
        module load samtools/1.10.1
        module load bedtools/2.29.2
        bwa mem -t {threads} {params.genome} {input.cleanFastQ1} {input.cleanFastQ2} > {params.sam}
        samtools view -Subh -o {params.bam} {params.sam}
        samtools sort -o {output.sortedBam} {params.bam}
        """

rule SS_index_bam:
    input:
        sortedBam = "../results/{dataset}/assembly/megahit_g1000/mapped_reads/{sample}.mapped.sorted.bam"
    output:
        bamIndex = "../results/{dataset}/assembly/megahit_g1000/mapped_reads/{sample}.mapped.sorted.bam.bai"
    shell:
        """
        module load samtools/1.10.1
        samtools index {input.sortedBam}
        """


rule SS_metabat2_depth:
    input:
        sortedBam = "../results/{dataset}/assembly/megahit_g1000/mapped_reads/{sample}.mapped.sorted.bam",
        contigs = "../results/{dataset}/assembly/megahit_g1000/{sample}.megahit_g1000.fa",
        bamIndex = "../results/{dataset}/assembly/megahit_g1000/mapped_reads/{sample}.mapped.sorted.bam.bai"
    output:
        depth_fi = "../results/{dataset}/assembly/megahit_g1000/metabat2/{sample}/{sample}_depth.txt"
    conda:
        "../envs/metabat2.yml"
    threads: 20
    shell:
        """
        jgi_summarize_bam_contig_depths \
        --outputDepth {output.depth_fi} \
        --percentIdentity 97 \
        --minContigLength 1000 \
        --minContigDepth 1.0 \
        --referenceFasta {input.contigs} {input.sortedBam}
        """

rule SS_metabat2_bin:
    input:
        contigs = "../results/{dataset}/assembly/megahit_g1000/{sample}.megahit_g1000.fa",
        depth_fi = "../results/{dataset}/assembly/megahit_g1000/metabat2/{sample}/{sample}_depth.txt"
    output:
        bin_one = "../results/{dataset}/assembly/megahit_g1000/metabat2/{sample}/bins/bin.1.fa",
        bin_dir = directory("../results/{dataset}/assembly/megahit_g1000/metabat2/{sample}/bins")
    conda:
        "../envs/metabat2.yml"
    threads: 20
    shell:
        """
        metabat2 -t {threads} \
        --inFile {input.contigs} \
        --outFile {output.bin_dir}/bin \
        --abdFile {input.depth_fi}
        """
    

rule SS_checkm_analysis:
    input:
        bin_dir = directory("../results/{dataset}/assembly/megahit_g1000/metabat2/{sample}/bins")
    output:
        checkm_dir = directory("../results/{dataset}/assembly/megahit_g1000/metabat2/checkm/{sample}"),
        checkm_fi = "../results/{dataset}/assembly/megahit_g1000/metabat2/checkm/{sample}/{sample}_checkm_output.txt"
    threads: 12
    shell:
        """
        module load checkm/1.0.7
        checkm lineage_wf \
        --threads {threads} \
        --extension 'fa' \
        --file {output.checkm_fi} \
        {input.bin_dir} {output.checkm_dir}
        """

# to do: add config data for checkm
# add checkm plots

rule SS_checkm_plots:
    input:
        checkm_dir = directory("../results/{dataset}/assembly/megahit_g1000/metabat2/checkm/{sample}"),
        bin_dir = directory("../results/{dataset}/assembly/megahit_g1000/metabat2/{sample}/bins")
    output:
        plots_dir = directory("../results/{dataset}/assembly/megahit_g1000/metabat2/checkm/{sample}/plots"),
        qa_plot = "../results/{dataset}/assembly/megahit_g1000/metabat2/checkm/{sample}/plots/bin_qa_plot.png"
    shell:
        """
        module load checkm/1.0.7
        checkm bin_qa_plot \
        --extension 'fa' \
        {input.checkm_dir} {input.bin_dir} {output.plots_dir}
        """

# add single sample binning
# add mash 
