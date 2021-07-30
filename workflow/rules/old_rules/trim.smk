
#### Functions to id local read sets, return str ####

def get_fqs(wildcards):
    return samples.loc[wildcards.sample, ["f_fq1", "f_fq2", "s_fq1", "s_fq2"]].dropna()

def get_fq_f_r1(wildcards):
    tmp = get_fqs(wildcards)
    return "../data/reads/" + tmp["f_fq1"]

def get_fq_f_r2(wildcards):
    tmp = get_fqs(wildcards)
    return "../data/reads/" + tmp["f_fq2"]

def get_fq_s_r1(wildcards):
    tmp = get_fqs(wildcards)
    return "../data/reads/" + tmp["s_fq1"]

def get_fq_s_r2(wildcards):
    tmp = get_fqs(wildcards)
    return "../data/reads/" + tmp["s_fq2"]
# Function to get sample id for reads, returns str
def get_sample(wildcards):
    return "../data/reads/" + samples.loc[wildcards.sample, ["sample", "f_fq1", "f_fq2", "s_fq1", "s_fq2"]].dropna()["sample"]


#### Rules to trim first and second round of sequencing data ###

rule first_fastp:
    input:
        f_r1=get_fq_f_r1,
        f_r2=get_fq_f_r2
    output:
        tf_r1=temp("../data/trimmed/first/{sample}.f.r1.fastq.gz"),
        tf_r2=temp("../data/trimmed/first/{sample}.f.r2.fastq.gz"),
        json="../data/trimmed/first/{sample}_fastp.json",
        html="../data/trimmed/first/{sample}_fastp.html"
    conda:
        "../envs/seq_processing.yml"
    threads: 16
    shell:
        "fastp -i {input.f_r1} -I {input.f_r2} --out1 {output.tf_r1} --out2 {output.tf_r2} --detect_adapter_for_pe --thread {threads} --length_required 50 -j {output.json} -h {output.html} -V"


rule second_fastp:
    input:
        s_r1=get_fq_s_r1,
        s_r2=get_fq_s_r2
    output:
        ts_r1=temp("../data/trimmed/second/{sample}.s.r1.fastq.gz"),
        ts_r2=temp("../data/trimmed/second/{sample}.s.r2.fastq.gz"),
        json="../data/trimmed/second/{sample}_fastp.json",
        html="../data/trimmed/second/{sample}_fastp.html"
    conda:
        "../envs/seq_processing.yml"
    threads: 16
    shell:
        "fastp -i {input.s_r1} -I {input.s_r2} --out1 {output.ts_r1} --out2 {output.ts_r2} --detect_adapter_for_pe --thread {threads} --length_required 50 -j {output.json} -h {output.html} -V"

