from snakemake.utils import validate
import pandas as pd

###########################################################
###    Takes pipeline configuration to set rule all     ###
###########################################################

def get_rules(wildcards):
    all_rules = []
    if config["FASTQC"]:
        all_rules.append("results/fastqc_out/multiqc_report.html")
    if config["MODULE_READ_QC"]:
        if config["MERGE_READS"]:
            all_rules = all_rules + expand(
                "results/bbmerge_out/{sample}/{sample}.bbmerge.fastq.gz",
                sample=samples["sample"])
        elif config["DECONVOLUTE"]:
            if config["BOWTIE2"]:
                    all_rules.append("results/bowtie_out/flagstat_summary.txt")
        elif config["COMPLEXITY_FILTER"]:
            all_rules = all_rules + expand(
                "results/bbduk_out/{sample}/{sample}.bbduk.r1.fastq.gz", 
                sample=samples["sample"])
            all_rules = all_rules + expand(
                "results/bbduk_out/{sample}/{sample}.bbduk.r2.fastq.gz", 
                sample=samples["sample"])
        elif config["TRIM_READS"]:
            all_rules = all_rules + expand(
                    "results/fastp_out/{sample}/{sample}.fastp.r1.fastq.gz", 
                    sample=samples["sample"])
            all_rules = all_rules + expand(
                    "results/fastp_out/{sample}/{sample}.fastp.r2.fastq.gz", 
                    sample=samples["sample"])
        else:
            pass
    if config["NONPAREIL"]:
        if config["MERGE_READS"]:
            all_rules = all_rules + expand("results/nonpareil.merged_out/{sample}/{sample}.npo",
                sample=samples["sample"])
        #else:
        all_rules = all_rules + expand(
            "results/nonpareil_out/{sample}/{sample}.npo", 
            sample=samples["sample"])
    if config["METAPHLAN"]:
        all_rules.append("results/metaphlan_bowtie_out/merged_metaphlan_profile_genus.tsv")
        all_rules.append("results/metaphlan_bowtie_out/merged_metaphlan_profile_species.tsv")
    if config["KRAKEN2"]:
        all_rules.append("results/kraken/merged_kraken_report_profile.tsv")
        #all_rules = all_rules + expand("results/kraken/{sample}/{sample}_kraken2out.txt", sample=samples["sample"])
    if config["METAXA2"]:
        all_rules = all_rules + expand("results/metaxa2/{sample}/{sample}_metaxa2.taxonomy.txt", sample=samples["sample"])
    if config["ASSEMBLE"]:
        #if config["MEGAHIT"]:
        #    all_rules = all_rules + expand(
        #        "results/megahit_out/{sample}/{sample}.contigs.fa", 
        #        sample=samples["sample"])
        #if config["SPADES"]:
        #    all_rules = all_rules + expand(
        #        "results/spades_out/{sample}/scaffolds.fasta", 
        #        sample=samples["sample"])
        if config["SPADES"] or config["MEGAHIT"]:
            all_rules.append("results/quast_out/multiqc_report.html")
            #all_rules = all_rules + expand(
            #    "results/{assembler}_parsed/{sample}/{sample}.parsed_assembly.fa", 
            #    sample=samples["sample"],
            #    assembler=ASSEMBLER)
    if config["METABAT2"]:
        metabat2_results = directory(expand(
            "results/metabat_{assembler}_out/{sample}/bins/", 
            assembler=ASSEMBLER, #["megahit", "spades"]
            sample=samples["sample"]))
        all_rules = all_rules + metabat2_results
    return all_rules

###########################################################
###  Helper functions for determining readQC rule I/O   ###
###########################################################

### Helper functions for getting initial reads ###

def get_read_path_v2(wildcards):
    return samples.loc[wildcards.sample, ["sample", "dataset", "forward_read", "reverse_read"]].dropna()

def get_r1(wildcards):
    tmp = get_read_path_v2(wildcards)
    return tmp["forward_read"]

def get_r2(wildcards):
    tmp = get_read_path_v2(wildcards)
    return tmp["reverse_read"]

### Helper functions for getting trimmed reads ###

def get_trimmed_read1(wildcards):
    return "results/fastp_out/{sample}/{sample}.fastp.r1.fastq.gz"

def get_trimmed_read2(wildcards):
    return "results/fastp_out/{sample}/{sample}.fastp.r2.fastq.gz"

### Helper functions for getting complexity-filtered reads ###

def get_complex_read1(wildcards):
    return "results/bbduk_out/{sample}/{sample}.bbduk.r1.fastq.gz"

def get_complex_read2(wildcards):
    return "results/bbduk_out/{sample}/{sample}.bbduk.r2.fastq.gz"

### Helper functions for getting control-decontaminated reads ###

def get_control_decontaminated_read1(wildcards):
    return "results/negative_out/{sample}/{sample}.clean.r1.fastq.gz"

def get_control_decontaminated_read2(wildcards):
    return "results/negative_out/{sample}/{sample}.clean.r2.fastq.gz"

### Helper functions for getting host-decontaminated reads ###

def get_decontaminated_read1(wildcards):
    return "results/bowtie_out/{sample}/{sample}.bowtie.r1.fastq.gz"

def get_decontaminated_read2(wildcards):
    return "results/bowtie_out/{sample}/{sample}.bowtie.r2.fastq.gz"

### Helper functions for getting merged reads ###

def get_merged_reads(wildcards):
    return  "results/bbmerge_out/{sample}/{sample}.bbmerge.fastq.gz"

### Helper functions for determining inputs for control decontamination module ###

def get_input_control_decontaminated_read1(wildcards):
    if config["COMPLEXITY_FILTER"]:
        return get_complex_read1(wildcards)
    elif config["TRIM_READS"]:
        return get_trimmed_read1(wildcards)
    else:
        return get_r1(wildcards)

def get_input_control_decontaminated_read2(wildcards):
    if config["COMPLEXITY_FILTER"]:
        return get_complex_read2(wildcards)
    elif config["TRIM_READS"]:
        return get_trimmed_read2(wildcards)
    else:
        return get_r2(wildcards)

### Helper functions for determining penultimate paired/merged reads ###

def get_penultimate_read1(wildcards):
    if config["CTRL_CONTAMINANT_FILTER"]:
        return get_control_decontaminated_read1(wildcards)
    if config["COMPLEXITY_FILTER"]:
        return get_complex_read1(wildcards)
    elif config["TRIM_READS"]:
        return get_trimmed_read1(wildcards)
    else:
        return get_r1(wildcards)

def get_penultimate_read2(wildcards):
    if config["CTRL_CONTAMINANT_FILTER"]:
        return get_control_decontaminated_read2(wildcards)
    if config["COMPLEXITY_FILTER"]:
        return get_complex_read2(wildcards)
    elif config["TRIM_READS"]:
        return get_trimmed_read2(wildcards)
    else:
        return get_r2(wildcards)
### Helper functions for determining final paired/merged reads ###

def get_final_read1(wildcards):
    if config["DECONVOLUTE"]:
        return get_decontaminated_read2(wildcards)
    elif config["CTRL_CONTAMINANT_FILTER"]:
        return get_control_decontaminated_read1(wildcards)
    elif config["COMPLEXITY_FILTER"]:
        return get_complex_read1(wildcards)
    elif config["TRIM_READS"]:
        return get_trimmed_read1(wildcards)
    else:
        return get_r1(wildcards)

def get_final_read2(wildcards):
    if config["DECONVOLUTE"]:
        return get_decontaminated_read2(wildcards)
    elif config["CTRL_CONTAMINANT_FILTER"]:
        return get_control_decontaminated_read2(wildcards)
    elif config["COMPLEXITY_FILTER"]:
        return get_complex_read2(wildcards)
    elif config["TRIM_READS"]:
        return get_trimmed_read2(wildcards)
    else:
        return get_r2(wildcards)

def get_final_merged_read(wildcards):
    return get_merged_reads(wildcards)

###########################################################
###   Helper functions for determining misc. rule I/O   ###
###########################################################

### Helper functions for configuring quast multiqc input ###
def get_multiqc_input():
    """
    Gets list of files that were expected to be made based on config settings and 
    creates a list of fastqc output files to expect for multiqc input
    """
    multiqc_in = []
    if config["MERGE_READS"]:
        multiqc_in = multiqc_in + expand("results/fastqc_out/bbmerge_qc/{sample}/{sample}.bbmerge_fastqc.html",
            sample=samples["sample"])
    if config["CTRL_CONTAMINANT_FILTER"]:
        multiqc_in = multiqc_in + expand("results/fastqc_out/bowtie_qc/{sample}/{sample}.clean.r1_fastqc.html",
            sample=samples["sample"])
    if config["DECONVOLUTE"]:
        multiqc_in = multiqc_in + expand("results/fastqc_out/bowtie_qc/{sample}/{sample}.bowtie.r1_fastqc.html",
            sample=samples["sample"])
    if config["COMPLEXITY_FILTER"]:
        multiqc_in = multiqc_in + expand("results/fastqc_out/bbduk_qc/{sample}/{sample}.bbduk.r1_fastqc.html",
            sample=samples["sample"])
    if config["TRIM_READS"]:
        multiqc_in = multiqc_in + expand("results/fastqc_out/fastp_qc/{sample}/{sample}.fastp.r1_fastqc.html",
            sample=samples["sample"])
    multiqc_in = multiqc_in + expand("results/fastqc_out/raw_qc/{sample}/{sample}.raw.r1_fastqc.html",
        sample=samples["sample"])
    return multiqc_in

### Helper functions for configuring quast multiqc input ###

def get_multiqc_quast_input(wildcards):
    #if config["MERGE_READS"] and config["SPADES"]:
    #    return expand("results/quast_out/spades.merged/{sample}/report.html", 
    #        sample=samples["sample"],
    #        assembler=ASSEMBLER)
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
controls = config["controls"]
negative_controls = config["negative_controls"]

ASSEMBLER = get_assemblers()
