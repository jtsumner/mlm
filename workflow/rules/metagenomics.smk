import glob
import pandas as pd
from snakemake.utils import validate


### Run FastP to trim/filter reads ###

rule fastp_pe:
    input:
        r1 = get_r1,
        r2 = get_r2
    output:
        r1Filtered = "../results/{dataset}/filtered/{sample}.filtered.R1.fastq.gz",
        r2Filtered = "../results/{dataset}/filtered/{sample}.filtered.R2.fastq.gz",
        json = "../results/{dataset}/filtered/{sample}_fastp.json",
        html = "../results/{dataset}/filtered/{sample}_fastp.html"
    conda:
        "../envs/seq_processing.yml"
    threads: 16
    shell: 
        "fastp -i {input.r1} -I {input.r2} --out1 {output.r1Filtered} --out2 {output.r2Filtered} --detect_adapter_for_pe --thread {threads} --length_required 50 -j {output.json} -h {output.html} -V"


rule fastqc:
    input: 
        "../results/{dataset}/filtered/{sample}.filtered.R1.fastq.gz",
        "../results/{dataset}/filtered/{sample}.filtered.R2.fastq.gz"
    output:
        "../results/{dataset}/filtered/fastqc/{sample}.filtered.R1_fastqc.html",
        "../results/{dataset}/filtered/fastqc/{sample}.filtered.R2_fastqc.html"
    params:
        outDir = "../results/{dataset}/filtered/fastqc"
    threads: 12
    shell:
        "module load fastqc/0.11.5 ; fastqc -t {threads} {input} --outdir {params.outDir}"


### Remove contaminant reads aligning to human reference genome ###

rule bwa_map:
    input:
        r1Filtered = "../results/{dataset}/filtered/{sample}.filtered.R1.fastq.gz",
        r2Filtered = "../results/{dataset}/filtered/{sample}.filtered.R2.fastq.gz"
    output:
        sam = temp("../results/{dataset}/bwa/{sample}.mapped.sam"),
        bam = "../results/{dataset}/bwa/{sample}.mapped.bam",
        sortedBam = "../results/{dataset}/bwa/{sample}.mapped.sorted.bam",
        unmappedBam = "../results/{dataset}/bwa/{sample}.unmapped.bam",
        cleanFastQ1 = "../results/{dataset}/bwa/{sample}.clean.R1.fastq",
        cleanFastQ2 = "../results/{dataset}/bwa/{sample}.clean.R2.fastq"
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


### Setup Metaphlan. Run Metaphlan on samples to make abundance tables ###

def metaphlan_merge_inputs(wildcards):
    files = expand("../results/{dataset}/abundance/metaphlan/{sample}.metaphlan_profile.txt",
        zip, sample=samples["sample"], dataset=samples["dataset"])
    return files


rule metaphlan_setup:
    output:
        metaphlan_db = directory("../resources/metaphlan_db")
    conda: 
        "../envs/metaphlan.yml"
    params:
        metaphlan_idx = config["metaphlan_idx"] # Index for metaphlan
    threads: 10
    shell:
        """
        metaphlan --install --index {params.metaphlan_idx} --bowtie2db {output.metaphlan_db} --nproc {threads}
        """


rule metaphlan:
    input:
        metaphlan_db = rules.metaphlan_setup.output.metaphlan_db,
        cleanFastQ1 = "../results/{dataset}/bwa/{sample}.clean.R1.fastq",
        cleanFastQ2 = "../results/{dataset}/bwa/{sample}.clean.R2.fastq"
    output:
        profile = "../results/{dataset}/abundance/metaphlan/{sample}.metaphlan_profile.txt",
        bowtie_out = "../results/{dataset}/abundance/metaphlan/{sample}.bowtie2.bz2"
    conda: 
        "../envs/metaphlan.yml"
    params:
        metaphlan_idx = config["metaphlan_idx"] # Index for metaphlan
    threads: 20
    shell:
        """
        metaphlan {input.cleanFastQ1},{input.cleanFastQ2} \
        --bowtie2out {output.bowtie_out} \
        --index {params.metaphlan_idx} \
        --bowtie2db {input.metaphlan_db} \
        --nproc {threads} \
        --input_type fastq \
        -o {output.profile}
        """


rule metaphlan_merge:
    input:
        metaphlan_merge_inputs
    output:
        "../results/allDatasets/metaphlan/merged_abundance_table.allDatasets.txt"
    conda:
        "../envs/metaphlan.yml"
    shell:
        """
        merge_metaphlan_tables.py {input} > {output}
        """


rule metaphlan_species_abundance:
    input:
        "../results/allDatasets/metaphlan/merged_abundance_table.allDatasets.txt"
    output:
        "../results/allDatasets/metaphlan/merged_abundance_table.species.allDatasets.txt"
    conda:
        "../envs/metaphlan.yml"
    shell:
        """
        grep -E "s__|clade" {input} | sed 's/^.*s__//g' \
        | cut -f1,3- | sed -e 's/clade_name/sample/g' > {output}
        """


rule metaphlan_genus_abundance:
    input:
        "../results/allDatasets/metaphlan/merged_abundance_table.allDatasets.txt"
    output:
        "../results/allDatasets/metaphlan/merged_abundance_table.genus.allDatasets.txt"
    conda:
        "../envs/metaphlan.yml"
    shell:
        """
        grep -E "g__|clade" {input} | sed 's/^.*g__//g' \
        | grep -v s__ |cut -f1,3- | sed -e 's/clade_name/sample/g' > {output}
        """


rule hclust:
    input:
        "../results/allDatasets/metaphlan/merged_abundance_table.species.allDatasets.txt"
    output:
        "../results/allDatasets/metaphlan/abundance_heatmap_species.allDatasets.png"
    conda:
        "../envs/hclust.yml"
    shell:
        """
        hclust2.py -i {input} -o {output} --f_dist_f braycurtis --s_dist_f braycurtis --cell_aspect_ratio 0.5 -l --flabel_size 10 --slabel_size 10 --max_flabel_len 100 --max_slabel_len 100 --minv 0.1 --dpi 300
        """

### Setup and Execute Kaiju ###
 
rule kaiju_setup:
    output:
        tar = "../resources/kaiju_head/kaiju-v1.8.0-linux-x86_64.tar.gz"
    threads: 10
    params:
        kaiju_archive = "https://github.com/bioinformatics-centre/kaiju/releases/download/v1.8.0/kaiju-v1.8.0-linux-x86_64.tar.gz",
        kaiju_old_dir = "kaiju-v1.8.0-linux-x86_64-static",
        kaiju_head = "../resources/kaiju_head",
        kaijuDB = "../resources/kaiju_head/kaijuDB",
        database = "refseq"
    shell:
        """
        wget {params.kaiju_archive} -P {params.kaiju_head}
        tar -xvzf {params.kaiju_head}/kaiju-v1.8.0-linux-x86_64.tar.gz -C {params.kaiju_head} --transform s/{params.kaiju_old_dir}/kaijuDir/
        cd {params.kaijuDB}
        {params.kaiju_head}/kaijuDir/kaiju-makedb -s {params.database} -t {threads}
        """
"""
rule kaiju_db:
    input:
    output:
    threads:
    shell:

rule kaiju_refseq:
    input:
    output:
    threads:
    shell:

"""


"""
rule metaphlan_abundance:
"""

"""
        # samtools view -Subh {output.sam} | samtools sort –o {output.bam} -

rule multiqc:
    input:
        expand("../results/{dataset}/filtered/{sample}.filtered.R1.fastq.gz"),
        expand("../results/{dataset}/filtered/{sample}.filtered.R2.fastq.gz")
    output:
        expand("../results/{dataset}/filtered/fastqc/multiqc_report.html")
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
