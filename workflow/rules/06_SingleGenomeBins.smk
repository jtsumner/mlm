
import glob
import pandas as pd
from snakemake.utils import validate

### For single sample assemblies 

rule index_contigs:
    input:
        "results/{assembler}_parsed/{sample}/{sample}.fa"
    output:
        "results/{assembler}_parsed/{sample}/{sample}.1.bt2"
    params:
        index_name = "results/{assembler}_parsed/{sample}/{sample}"
    shell:
        """
        module load bowtie2/2.4.5
        module load samtools/1.10.1

        bowtie2-build {input} {params.index_name}
        """
# TODO add reduced time to index contis ... 

rule map2contigs:
    input:
        r1_clean = get_final_read1,
        r2_clean = get_final_read2,
        parsed_contigs = "results/{assembler}_parsed/{sample}/{sample}.fa",
        index = "results/{assembler}_parsed/{sample}/{sample}.1.bt2"
    output:
        sorted_bam = "results/{assembler}_parsed/{sample}/{sample}.sorted.bam",
        flagstat = "results/{assembler}_parsed/{sample}/{sample}.flagstat.tsv",
        bam_index = "results/{assembler}_parsed/{sample}/{sample}.sorted.bam.bai"
    threads: 1
    params:
        index_name = "results/{assembler}_parsed/{sample}/{sample}"
    benchmark:
        "results/logs/benchmarks/{assembler}_{sample}.map2contigs_benchmark.txt"
    threads: 10
    resources:
        mem="15G",
        time="01:00:00"
    shell:
        """
        module load bowtie2/2.4.5
        module load samtools/1.10.1

        bowtie2 -p {threads} -x {params.index_name} --very-sensitive-local -1 {input.r1_clean} -2 {input.r2_clean}| \
        samtools view -bS -@ {threads}| \
        samtools sort -@ {threads} -n -o {output.sorted_bam}

        samtools index {output.sorted_bam}

        samtools flagstat -@ {threads} -O tsv {output.sorted_bam} > {output.flagstat}
        """

# TODO take out sort in index bam and remove -n in map2contigs + put bams spades_parse
# rule index_bam:
#     input:
#         bam_sorted = "results/{assembler}_parsed/{sample}/{sample}.sorted.bam"
#     output:
#         bam_index = "results/{assembler}_parsed/{sample}/{sample}.sorted.bam.bai"
#     threads: 1
#     resources:
#         time = "00:05:00"
#     shell:
#         """
#         module load samtools/1.10.1
#         samtools sort -@ {threads} -o {input.bam_sorted} {input.bam_sorted}
#         samtools index {input.bam_sorted}
#         """

rule metabat_depth:
    input:
        sorted_bam = "results/{assembler}_parsed/{sample}/{sample}.sorted.bam",
        contigs = "results/{assembler}_parsed/{sample}/{sample}.fa",
        bam_index = "results/{assembler}_parsed/{sample}/{sample}.sorted.bam.bai"
    output:
        depth_fi = "results/metabat_{assembler}_out/{sample}/{sample}_depth.txt"
    conda:
        "../envs/metabat2.yml"
    threads: 1
    resources:
        time = "00:05:00"
    shell:
        """
        jgi_summarize_bam_contig_depths \
            --outputDepth {output.depth_fi} \
            --referenceFasta {input.contigs} {input.sorted_bam}
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
    threads: 2
    resources:
        mem="10g",
        time="00:05:00"
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
