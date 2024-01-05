import glob
import pandas as pd
from snakemake.utils import validate

############################
###   PART 1A: MEGAHIT   ###
############################

rule megahit:
    input:
        r1_clean = get_final_read1,
        r2_clean = get_final_read2
    output:
        scaffolds = "results/megahit_out/{sample}/{sample}.contigs.fa",
    params:
        out_dir = "results/megahit_out/{sample}"
    conda:
        "../envs/megahit.yml"
    threads: 20
    resources:
        mem="100g",
        time="00:30:00",
        partition="genomics",
        account="b1042"
    shell:
        """
        megahit -t {threads} \
            -m 0.9 -1 {input.r1_clean} \
            -2 {input.r2_clean} \
            --out-prefix {wildcards.sample} \
            -o {params.out_dir}_tmp
        mv {params.out_dir}_tmp/* {params.out_dir}
        rmdir {params.out_dir}_tmp
        """
        
############################
###   PART 1B: SPADES    ###
############################

rule spades:
    input:
        r1_clean = get_final_read1,
        r2_clean = get_final_read2
    output:
        scaffolds="results/spades_out/{sample}/scaffolds.fasta"
    params:
        out_dir=directory("results/spades_out/{sample}")

    threads: 20
    resources:
        mem="200g",
        time="03:00:00",
        partition="genomics",
        account="b1042"
    shell:
        """
        module load spades/3.14.1
        spades.py -1 {input.r1_clean} -2 {input.r2_clean} -o {params.out_dir} -t {threads} -m 300 --meta
        """
# -k 21,33,55,77,99,127 --only-assembler

############################
##  PART 1C: COASSEMBLY   ##
############################


def get_final_reads_all1(wildcards):
    r1 = get_final_read1(wildcards)
    return expand(r1, sample=samples[~samples["sample"].isin(controls)]["sample"])

def get_final_reads_all2(wildcards):
    r2 = get_final_read2(wildcards)
    return expand(r2, sample=samples[~samples["sample"].isin(controls)]["sample"])

rule concat_reads:
    input:
        r1_clean = get_final_reads_all1,
        r2_clean = get_final_reads_all2
    output:
        r1_concat = "results/coassembly_out/concat_reads.clean.r1.fastq.gz",
        r2_concat = "results/coassembly_out/concat_reads.clean.r2.fastq.gz"
    shell:
        """
        cat {input.r1_clean} > {output.r1_concat}
        cat {input.r2_clean} > {output.r2_concat}
        """

# use rule megahit as megahit_coassembly with:
#     input:
#         r1_concat = "results/coassembly_out/concat_reads.clean.r1.fastq.gz",
#         r2_concat = "results/coassembly_out/concat_reads.clean.r1.fastq.gz"
#     output:
#         scaffolds = "results/coassembly_out/{sample}.contigs.fa"
#     params:
#         out_dir = "results/coassembly_out"
#     threads: 60
#     resources:
#         mem="200g",
#         time="10:00:00"



############################
#   PART 2A: QC FILTER     #
############################

rule drop_short_contigs_spades:
    input:
        scaffolds = "results/spades_out/{sample}/scaffolds.fasta"
    output:
        parsed = "results/spades_parsed/{sample}/{sample}.fa"
    params:
        min_length = "200",
        assembler="spades"
    
    conda:
        "../envs/biopython.yml"
    shell:
        """
        .snakemake/conda/3cac06327597fe89b28a87b1f6e65a77_/bin/python workflow/scripts/parse_contigs.py --sample {wildcards.sample} \
            --scaffolds {input.scaffolds} \
            --parsed {output.parsed} \
            --min_length {params.min_length} \
            --assembler {params.assembler}
        """

use rule drop_short_contigs_spades as drop_short_contigs_megahit with:
    input:
        scaffolds = "results/megahit_out/{sample}/{sample}.contigs.fa"
    output:
        parsed = "results/megahit_parsed/{sample}/{sample}.fa"
    params:
        min_length = "200",
        assembler="megahit"

############################
#    PART 2B: QC QUAST     #
############################

rule quast_megahit:
    input:
        scaffolds = "results/megahit_out/{sample}/{sample}.contigs.fa"
    output:
        out_dir=directory("results/quast_out/megahit/{sample}"),
        report="results/quast_out/megahit/{sample}/report.html"
    params:
        out_dir="results/quast_out/megahit/{sample}"
    threads: 1
    conda:
        "../envs/genome_qc.yml"
    shell:
        """
        quast.py -o {params.out_dir} --threads {threads} --min-contig 500 -L {input}
        """

rule quast_spades:
    input:
        scaffolds="results/spades_out/{sample}/scaffolds.fasta"
    output:
        report="results/quast_out/spades/{sample}/report.html"
    params:
        out_dir="results/quast_out/spades/{sample}"
    threads: 1
    conda:
        "../envs/genome_qc.yml"
    shell:
        """
        quast.py -o {params.out_dir} --threads {threads} --min-contig 500 -L {input}
        """

rule multiqc_quast:
    input:
        quast_reports = get_multiqc_quast_input
    output:
        multiqc_report = "results/quast_out/multiqc_report.html"
    params:
        out_dir="results/quast_out"
    shell:
        """
        module load multiqc
        multiqc --outdir {params.out_dir} --dirs --dirs-depth 2 results/quast_out/ -f
        """

############################
#   PART 3A: ANNOTATION    #
############################

rule prep_annotation:
    """
    Takes spades headers and transforms them to be >NODE_#### only
    """
    input:
        scaffolds="results/spades_out/{sample}/scaffolds.fasta"
    output:
        scaffolds_clean="results/prokka_out/tmp_scaffolds/{sample}_scaffolds_clean.fasta"
    threads: 1
    resources:
        mem="3g",
        time="00:05:00"
    shell:
        """
        awk '/^>/ {{ sub(/_length_[0-9]+_cov_[0-9.]+/, "", $0) }} 1' {input.scaffolds} > {output.scaffolds_clean}
        """

rule annotate_prokka:
    input:
        scaffolds_clean="results/prokka_out/tmp_scaffolds/{sample}_scaffolds_clean.fasta"
    output:
        annotation="results/prokka_out/{sample}/{sample}.tsv",
        amino_acids="results/prokka_out/{sample}/{sample}.faa",
        out_dir=directory("results/prokka_out/{sample}/")
    threads: 20
    resources:
        mem="25g",
        time="01:00:00"
    conda:
        "../envs/prokka.yml"
    shell:
        """
        module load prokka
        prokka {input.scaffolds_clean} --cpus {threads} \
            --metagenome \
            --locustag {wildcards.sample} \
            --prefix {wildcards.sample} \
            --outdir {output.out_dir}/ \
            --addgenes \
            --mincontiglen 200 \
            --force
        """

rule merge_prokka_aa:
    """
    Takes prokka outputs from SAMPLES (see config to ID controls)
    and concatenates to make a master file
    """
    input:
        annotation=expand("results/prokka_out/{sample}/{sample}.faa", sample=samples[~samples["sample"].isin(controls)]["sample"])
    output:
        merged_annotation="results/prokka_out/merged_fastas.faa"
    threads: 1
    resources:
        mem="3g",
        time="00:30:00"
    shell:
        """
        cat {input.annotation} > {output.merged_annotation}
        """

rule eggnog_setup:
    output:
        eggnog_db = "resources/eggnog_db/eggnog_db.db"
    conda:
        "../env/eggnog.yml"
    shell:
        """
        download_eggnog_data.py -y --data_dir resources/eggnog_db/
        """

rule emapper:
    input:
        merged_annotation="results/prokka_out/merged_fastas.faa",
        eggnog_db = "resources/eggnog_db/eggnog_db.db"
    output:
        annotations = "results/eggnog/merged_annotation.emapper.annotations"
    params:
        prefix = "merged_annotations"
    threads: 20
    resources:
        mem="30G"
    conda:
        "../envs/eggnog.yml"
    shell:
        """
        emapper.py \
            -m diamond \
            --cpu {threads} \
            -i {input.merged_annotation} \
            -o  {params.prefix} \
            --out_dir $(dirname {output.annotations}) \
            --data_dir $(dirname {input.eggnog.db})
        """
