
import glob
import pandas as pd
from snakemake.utils import validate

### For single sample assemblies 

rule index_sample_contigs:
    input:
        "../results/allDatasets/single_sample_assemblies/allSamples.megahit_g1000.fa"
    output:
        "../results/allDatasets/single_sample_assemblies/megahit_g1000.bwt"
    shell:
        """
        module load bwa/0.7.17
        module load samtools/1.10.1

        bwa index -p ../results/allDatasets/single_sample_assemblies/megahit_g1000 {input}
        """


rule map_reads_to_assembly:
    input:
        cleanFastQ1 = "../results/{dataset}/bwa/{sample}.clean.R1.fastq",
        cleanFastQ2 = "../results/{dataset}/bwa/{sample}.clean.R2.fastq",
        index = "../results/allDatasets/single_sample_assemblies/megahit_g1000.bwt"
    output:
        sortedBam = "../results/allDatasets/single_sample_assemblies/mapped_reads/{dataset}/{sample}.mapped.sorted.bam"
    params:
        genome = "../results/allDatasets/single_sample_assemblies/megahit_g1000",
        sam = "../results/allDatasets/single_sample_assemblies/mapped_reads/{dataset}/{sample}.mapped.sam",
        bam = "../results/allDatasets/single_sample_assemblies/mapped_reads/{dataset}/{sample}.mapped.bam",
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

rule index_bam:
    input:
        sortedBam = "../results/allDatasets/single_sample_assemblies/mapped_reads/{dataset}/{sample}.mapped.sorted.bam"
    output:
        bamIndex = "../results/allDatasets/single_sample_assemblies/mapped_reads/{dataset}/{sample}.mapped.sorted.bam.bai"
    shell:
        """
        module load samtools/1.10.1
        samtools index {input.sortedBam}
        """


rule metabat2_depth:
    input:
        sortedBam = expand("../results/allDatasets/single_sample_assemblies/mapped_reads/{dataset}/{sample}.mapped.sorted.bam",
            zip, sample=samples["sample"], dataset=samples["dataset"]),
        contigs = "../results/allDatasets/single_sample_assemblies/allSamples.megahit_g1000.fa",
        bamIndex = expand("../results/allDatasets/single_sample_assemblies/mapped_reads/{dataset}/{sample}.mapped.sorted.bam.bai",
            zip, sample=samples["sample"], dataset=samples["dataset"])
    output:
        depth_fi = "../results/allDatasets/single_sample_assemblies/metabat2/allSamples.megahit_g1000.fa.depth.txt"
    shell:
        """
        jgi_summarize_bam_contig_depths \
        --outputDepth {output.depth_fi} \
        --percentIdentity 97 \
        --minContigLength 1000 \
        --minContigDepth 1.0 \
        --referenceFasta {input.contigs} {input.sortedBam}
        """

rule metabat2_bin:
    input:
        contigs = "../results/allDatasets/single_sample_assemblies/allSamples.megahit_g1000.fa",
        depth_fi = "../results/allDatasets/single_sample_assemblies/metabat2/allSamples.megahit_g1000.fa.depth.txt"
    output:
        bin_one = "../results/allDatasets/single_sample_assemblies/metabat2/bins/bin.1.fa",
        bin_dir = directory("../results/allDatasets/single_sample_assemblies/metabat2/bins/")
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
    
rule bin_contigs:
    input:
        sortedBam = expand("../results/allDatasets/single_sample_assemblies/mapped_reads/{dataset}/{sample}.mapped.sorted.bam",
            zip, sample=samples["sample"], dataset=samples["dataset"]),
        contigs = "../results/allDatasets/single_sample_assemblies/allSamples.megahit_g1000.fa",
        bamIndex = expand("../results/allDatasets/single_sample_assemblies/mapped_reads/{dataset}/{sample}.mapped.sorted.bam.bai",
            zip, sample=samples["sample"], dataset=samples["dataset"])
    output:
        megahit_fi = "allSamples.megahit_g1000.fa.depth.txt"
    conda:
        "../envs/metabat2.yml"
    threads: 20
    shell:
        """
        runMetaBat.sh -t {threads} {input.contigs} {input.sortedBam} 
        """
    

rule checkm:
    input:
        inDir = directory("../results/allDatasets/single_sample_assemblies/megahit_genomeBins")
    output:
        outDir = directory("../results/allDatasets/single_sample_assemblies/megahit_genomeBins/checkm")
    threads: 12
    shell:
        """
        module load checkm/1.0.7
        checkm lineage_wf {input.inDir} -o 
        """

