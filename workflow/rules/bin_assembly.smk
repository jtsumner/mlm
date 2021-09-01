
import glob
import pandas as pd
from snakemake.utils import validate


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
        cleanFastQ2 = "../results/{dataset}/bwa/{sample}.clean.R2.fastq"
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

