
import glob
import pandas as pd
from snakemake.utils import validate

### For single sample assemblies 

rule index_contigs:
    input:
        "results/{assembler}_parsed/{sample}/{sample}.fa"
    output:
        "results/{assembler}_parsed/{sample}/{sample}.fa.bwt"
    shell:
        """
        module load bwa/0.7.17
        module load samtools/1.10.1

        bwa index {input}
        """


rule map2contigs:
    input:
        r1_clean = get_final_read1,
        r2_clean = get_final_read2,
        parsed_contigs = "results/{assembler}_parsed/{sample}/{sample}.fa",
        index = "results/{assembler}_parsed/{sample}/{sample}.fa.bwt"
    output:
        bam_sorted = "results/{assembler}_bams/{sample}/{sample}.sorted.bam"
    params:
        sam = "results/{assembler}_bams/{sample}/{sample}.sam",
        bam_unsorted = "results/{assembler}_bams/{sample}/{sample}.bam"
    benchmark:
        "logs/benchmarks/{assembler}_{sample}.map2contigs_benchmark.txt"
    threads: 15
    resources:
        mem="40G"
    shell:
        """
        module purge all
        module load bwa/0.7.17
        module load samtools/1.10.1
        module load bedtools/2.29.2
        bwa mem -t {threads} {input.parsed_contigs} {input.r1_clean} {input.r2_clean} > {params.sam}
        samtools view -Subh -o {params.bam_unsorted} {params.sam}
        samtools sort -o {output.bam_sorted} {params.bam_unsorted}
        """

rule index_bam:
    input:
        bam_sorted = "results/{assembler}_bams/{sample}/{sample}.sorted.bam"
    output:
        bam_index = "results/{assembler}_bams/{sample}/{sample}.sorted.bam.bai"
    shell:
        """
        module load samtools/1.10.1

        samtools index {input.bam_sorted}
        """

rule metabat_depth:
    input:
        bam_sorted = "results/{assembler}_bams/{sample}/{sample}.sorted.bam",
        contigs = "results/{assembler}_parsed/{sample}/{sample}.fa",
        bam_index = "results/{assembler}_bams/{sample}/{sample}.sorted.bam.bai"
    output:
        depth_fi = "results/metabat_{assembler}_out/{sample}/{sample}_depth.txt"
    conda:
        "../envs/metabat2.yml"
    threads: 10
    shell:
        """
        jgi_summarize_bam_contig_depths \
            --outputDepth {output.depth_fi} \
            --referenceFasta {input.contigs} {input.bam_sorted}
        """

rule metabat_bin:
    input:
        contigs = "results/{assembler}_parsed/{sample}/{sample}.fa",
        depth_fi = "results/metabat_{assembler}_out/{sample}/{sample}_depth.txt"
    output:
        bin_dir = directory("results/metabat_{assembler}_out/{sample}/bins")
        #bin_one = "results/metabat_{assembler}_out/{sample}/bins/bin.1.fa"
    conda:
        "../envs/metabat2.yml"
    threads: 5
    resources:
        mem="10g",
        time="00:30:00"
    shell:
        """
        metabat2 -t {threads} \
            --abdFile {input.depth_fi} \
            --inFile {input.contigs} \
            --outFile {output.bin_dir}/{wildcards.sample}.bin \
            --minContig 1500 \
            --seed=123456 \
            --unbinned \
            --verbose
        """

rule checkm_analysis:
    input:
        bin_dir = directory("results/metabat_{assembler}_out/{sample}/bins") #was spades
    output:
        checkm_dir = directory("results/checkm_{assembler}_out/{sample}"),
        checkm_fi = "results/checkm_{assembler}_out/{sample}/{sample}_checkm_output.txt"
    threads: 10
    resources:
        mem="80g",
        time="00:30:00"
    shell:
        """
        module load checkm/1.0.7
        checkm lineage_wf \
        --threads {threads} \
        --extension 'fa' \
        --file {output.checkm_fi} \
        --tab_table \
        {input.bin_dir} {output.checkm_dir}
        """

#use rule checkm_analysis as checkm_analysis_megahit with:
#    input:
#        bin_dir = directory("results/metabat_megahit_out/{sample}/bins")
#    output:
#        checkm_dir = directory("results/checkm_megahit_out/{sample}"),
#        checkm_fi = "results/checkm_megahit_out/{sample}/{sample}_checkm_output.txt"


# to do: add config data for checkm
# add checkm plots

rule SS_checkm_plots:
    input:
        checkm_dir = directory("results/metabat_out/checkm/{sample}"),
        bin_dir = directory("results/metabat_out/{sample}/bins")
    output:
        plots_dir = directory("results/metabat_out/checkm/{sample}/plots"),
        qa_plot = "results/metabat_out/checkm/{sample}/plots/bin_qa_plot.png"
    shell:
        """
        module load checkm/1.0.7
        checkm bin_qa_plot \
        --extension 'fa' \
        {input.checkm_dir} {input.bin_dir} {output.plots_dir}
        """

# add single sample binning
# add mash 
