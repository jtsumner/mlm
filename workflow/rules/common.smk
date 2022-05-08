from snakemake.utils import validate
import pandas as pd

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity

samples = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
samples.index.names = ["sample"]


def get_rules(wildcards):
    all_rules = []
    if config["FASTQC"]:
        all_rules = all_rules = all_rules + directory(
            expand(
                "results/fastqc_out/raw_qc/{sample}", 
                sample=samples["sample"] 
                )
        )

        if config["TRIM_READS"]:
            all_rules = all_rules = all_rules + expand(
                "results/fastqc_out/fastp_qc/{sample}/{sample}.fastp.r1_fastqc.html", 
                sample=samples["sample"]
            )
            all_rules = all_rules = all_rules + expand(
                "results/fastqc_out/fastp_qc/{sample}/{sample}.fastp.r2_fastqc.html", 
                sample=samples["sample"]
            )

        if config["ASSEMBLE"]:
            all_rules = all_rules = all_rules + expand(
                "results/fastqc_out/bwa_qc/{sample}/{sample}.fastp_bwa.r1_fastqc.html", 
                sample=samples["sample"]
            )
            all_rules = all_rules = all_rules + expand(
                "results/fastqc_out/bwa_qc/{sample}/{sample}.fastp_bwa.r2_fastqc.html", 
                sample=samples["sample"]
            )

    if config["TRIM_READS"]:
        all_rules = all_rules + expand(
            "results/fastp_out/{sample}/{sample}.fastp.r1.fastq.gz", 
            sample=samples["sample"]
        )
        all_rules = all_rules + expand(
            "results/fastp_out/{sample}/{sample}.fastp.r2.fastq.gz", 
            sample=samples["sample"] 
        )

    if config["DECONVOLUTE"]:
        all_rules = all_rules + expand(
            "results/bwa_out/{sample}/{sample}.fastp_bwa.r1.fastq", 
            sample=samples["sample"]
        )
        all_rules = all_rules + expand(
            "results/bwa_out/{sample}/{sample}.fastp_bwa.r2.fastq", 
            sample=samples["sample"]
        )

    if config["METAPHLAN"]:
        all_rules.append("results/metaphlan_merged/merged_metaphlan_hclust_species.png")
        all_rules.append("results/metaphlan_merged/merged_metaphlan_hclust_genus.png")
        #all_rules.append("results/metaphlan_merged/merged_metaphlan_unifrac_matrix.txt")

    if config["ASSEMBLE"]:
        all_rules = all_rules + expand(
            "results/megahit_out/{sample}/{sample}.contigs.fa", 
            sample=samples["sample"]
        )
        all_rules.append("results/quast_out/megahit/multiqc/multiqc_report.html")

        if config["SPADES"]:
            all_rules = all_rules + expand(
                "results/spades_out/{sample}/scaffolds.fasta", 
                sample=samples["sample"]
            )
            all_rules = all_rules + expand(
                "results/quast_out/spades/{sample}/report.html", 
                sample=samples["sample"]
            )
    if config["METABAT2"]:
        metabat2_results = expand(
            "results/{dataset}/assembly/megahit_g1000/metabat2/{sample}/bins/bin.1.fa", 
            sample=samples["sample"]
        )
        all_rules = all_rules + metabat2_results

    return all_rules


def get_read_path_v2(wildcards):
    return samples.loc[wildcards.sample, ["sample", "dataset", "forward_read", "reverse_read"]].dropna()


def get_r1(wildcards):
    tmp = get_read_path_v2(wildcards)
    return tmp["forward_read"]


def get_r2(wildcards):
    tmp = get_read_path_v2(wildcards)
    return tmp["reverse_read"]

