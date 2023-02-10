from snakemake.utils import validate
import pandas as pd


### Takes pipeline configuration to set rule all ###

def get_rules(wildcards):
    all_rules = []
    if config["FASTQC"]:
        if config["TRIM_READS"]:
            all_rules = all_rules = all_rules + expand(
                "results/fastqc_out/fastp_qc/{sample}/{sample}.fastp.r1_fastqc.html", 
                sample=samples["sample"]
            )
            all_rules = all_rules = all_rules + expand(
                "results/fastqc_out/fastp_qc/{sample}/{sample}.fastp.r2_fastqc.html", 
                sample=samples["sample"]
            )
        if config["DECONVOLUTE"]:
            pass
            #all_rules = all_rules = all_rules + expand(
            #    "results/fastqc_out/bwa_qc/{sample}/{sample}.fastp_bwa.r1_fastqc.html", 
            #    sample=samples["sample"])
            #all_rules = all_rules = all_rules + expand(
            #    "results/fastqc_out/bwa_qc/{sample}/{sample}.fastp_bwa.r2_fastqc.html", 
            #    sample=samples["sample"])
    if config["TRIM_READS"]:
        all_rules = all_rules + expand(
            "results/fastp_out/{sample}/{sample}.fastp.r1.fastq.gz", 
            sample=samples["sample"])
        all_rules = all_rules + expand(
            "results/fastp_out/{sample}/{sample}.fastp.r2.fastq.gz", 
            sample=samples["sample"])
    if config["DECONVOLUTE"]:
        if config["BWA"]:
            all_rules = all_rules + expand(
                "results/bwa_out/{sample}/{sample}.fastp_bwa.r1.fastq", 
                sample=samples["sample"])
            all_rules = all_rules + expand(
                "results/bwa_out/{sample}/{sample}.fastp_bwa.r2.fastq", 
                sample=samples["sample"])
        if config["BOWTIE2"]:
            all_rules = all_rules + expand(
                "results/bowtie_out/{sample}/{sample}.fastp_bowtie.r1.fastq.gz", 
                sample=samples["sample"])
            all_rules = all_rules + expand(
                "results/bowtie_out/{sample}/{sample}.fastp_bowtie.r2.fastq.gz", 
                sample=samples["sample"]) 
            all_rules.append("results/bowtie_out/flagstat_summary.txt")
        if config["NONPAREIL"]:
            all_rules = all_rules + expand(
                "results/nonpareil_out/{sample}/{sample}.npo", 
                sample=samples["sample"])
    if config["METAPHLAN"]:
        all_rules.append("results/metaphlan_bowtie_out/merged_metaphlan_profile_genus.tsv")
        all_rules.append("results/metaphlan_bowtie_out/merged_metaphlan_profile_species.tsv")
    if config["KRAKEN2"]:
        all_rules = all_rules + expand("results/kraken/{sample}/{sample}_kraken2out.txt", sample=samples["sample"])
    if config["METAXA2"]:
        all_rules = all_rules + expand("results/metaxa2/{sample}/{sample}_metaxa2.taxonomy.txt", sample=samples["sample"])
    if config["ASSEMBLE"]:
        if config["MEGAHIT"]:
            all_rules = all_rules + expand(
                "results/megahit_out/{sample}/{sample}.contigs.fa", 
                sample=samples["sample"])
        if config["SPADES"]:
            all_rules = all_rules + expand(
                "results/spades_out/{sample}/scaffolds.fasta", 
                sample=samples["sample"])
        if config["SPADES"] or config["MEGAHIT"]:
            all_rules.append("results/quast_out/multiqc_report.html")
            all_rules = all_rules + expand(
                "results/{assembler}_parsed/{sample}/{sample}.parsed_assembly.fa", 
                sample=samples["sample"],
                assembler=ASSEMBLER)
    if config["METABAT2"]:
        metabat2_results = directory(expand(
            "results/metabat_{assembler}_out/{sample}/bins/", 
            assembler=ASSEMBLER, #["megahit", "spades"]
            sample=samples["sample"]))
        all_rules = all_rules + metabat2_results
    return all_rules

### Helper functions for getting initial reads ###

def get_read_path_v2(wildcards):
    return samples.loc[wildcards.sample, ["sample", "dataset", "forward_read", "reverse_read"]].dropna()

def get_r1(wildcards):
    tmp = get_read_path_v2(wildcards)
    return tmp["forward_read"]

def get_r2(wildcards):
    tmp = get_read_path_v2(wildcards)
    return tmp["reverse_read"]

### Helper functions for determining assembly read input ###

def get_assembly_r1(wildcards):
    return "results/bowtie_out/{sample}/{sample}.fastp_bowtie.r1.fastq.gz"

def get_assembly_r2(wildcards):
    return "results/bowtie_out/{sample}/{sample}.fastp_bowtie.r2.fastq.gz"

### Helper functions for configuring quast multiqc input ###

def get_multiqc_quast_input(wildcards):
    if config["MEGAHIT"] and config["SPADES"]:
        quast_reports = expand(
            "results/quast_out/{assembler}/{sample}/report.html", 
            sample=samples["sample"],
            assembler=ASSEMBLER)
        return quast_reports
    elif config["MEGAHIT"] ^ config["SPADES"]:
        if config["MEGAHIT"]:
            quast_reports=expand(
                "results/quast_out/megahit/{sample}/report.html", 
                sample=samples["sample"])
        elif config["SPADES"]:
            quast_reports=expand(
                "results/quast_out/spades/{sample}/report.html", 
                sample=samples["sample"])
        return quast_reports
    else:
        print("CHECK that at least one assembler is set to True in config.yaml")

### Parse config file to set select assembler options as wildcard list ###

def get_assemblers():
    ASSEMBLER = []
    if config["MEGAHIT"]:
        ASSEMBLER.append("megahit")
    if config["SPADES"]:
        ASSEMBLER.append("spades")
    return ASSEMBLER

samples = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
samples.index.names = ["sample"]

ASSEMBLER = get_assemblers()
