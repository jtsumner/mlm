import glob
import pandas as pd
from snakemake.utils import validate


### Run FastP to trim/filter reads ###

# rule fastp_pe:
#     input:
#         r1 = get_r1,
#         r2 = get_r2
#     output:
#         r1Filtered = "results/fastp_out/{sample}/{sample}.fastp.r1.fastq.gz",
#         r2Filtered = "results/fastp_out/{sample}/{sample}.fastp.r2.fastq.gz",
#         json = "results/fastp_out/{sample}/{sample}_fastp.json",
#         html = "results/fastp_out/{sample}/{sample}_fastp.html"
#     conda:
#         "../envs/seq_processing.yml"
#     threads: 16
#     shell: 
#         "fastp -i {input.r1} -I {input.r2} --out1 {output.r1Filtered} --out2 {output.r2Filtered} --detect_adapter_for_pe --thread {threads} --length_required 50 -j {output.json} -h {output.html} -V"


# rule fastqc:
#     input: 
#         "results/fastp_out/{sample}/{sample}.fastp.r1.fastq.gz",
#         "results/fastp_out/{sample}/{sample}.fastp.r2.fastq.gz"
#     output:
#         "results/fastqc_out/fastp_qc/{sample}.fastp.r1_fastqc.html",
#         "results/fastqc_out/fastp_qc/{sample}.fastp.r2_fastqc.html"
#     params:
#         outDir = "results/fastqc_out/fastp_qc/"
#     threads: 12
#     shell:
#         "module load fastqc/0.11.5 ; fastqc -t {threads} {input} --outdir {params.outDir}"



# ### Remove contaminant reads aligning to human reference genome ###
# rule get_human_genome:
#     output:
#         "resources/genome/{}".format(config["genome_name"])
#     params:
#         human_genome = config["human_genome"]
#     shell:
#         """
#         cd resources/genome/
#         wget {params.human_genome}
#         """

# rule index_human_genome:
#     input:
#         "resources/genome/{}".format(config["genome_name"])
#     output:
#         "resources/genome/{}.ann".format(config["genome_name"])
#     params:
#         genome_index = config["genome_name"]
#     threads: 10
#     resources:
#         mem="20G",
#         time="04:00:00"
#     shell:
#         """
#         module load bwa/0.7.17
# 	    bwa index -a bwtsw {input} {params.genome_index}
#         """

# rule bwa_map:
#     input:
#         r1_filtered = "results/fastp_out/{sample}/{sample}.fastp.r1.fastq.gz",
#         r2_filtered = "results/fastp_out/{sample}/{sample}.fastp.r2.fastq.gz",
#         genome_idxed = "resources/genome/{}.ann".format(config["genome_name"]),
#         genome = "resources/genome/{}".format(config["genome_name"])
#     output:
#         r1_clean = "results/bwa_out/{sample}/{sample}.fastp_bwa.r1.fastq",
#         r2_clean = "results/bwa_out/{sample}/{sample}.fastp_bwa.r2.fastq"
#     params:
#         human_genome = config["genome_name"],
#         #genome = "resources/genome/{params.human_g}.ann",
#         sam = "results/bwa_out/{sample}/{sample}.mapped.sam",
#         bam = "results/bwa_out/{sample}/{sample}.mapped.bam",
#         sorted_bam = "results/bwa_out/{sample}/{sample}.mapped.sorted.bam",
#         unmapped_bam = "results/bwa_out/{sample}/{sample}.unmapped.bam"
#     threads: 20
#     resources:
#         mem="50G",
#         time="08:00:00"
#     shell:
#         """
#         module purge all
#         module load bwa/0.7.17
#         module load samtools/1.10.1
#         module load bedtools/2.29.2
#         bwa mem -t {threads} {input.genome} {input.r1_filtered} {input.r2_filtered} > {params.sam}
#         samtools view -Subh -o {params.bam} {params.sam}
#         samtools sort -o {params.sorted_bam} {params.bam}

#         samtools view -b -f 12 -F 256 -o {params.unmapped_bam} {params.sorted_bam}
#         bedtools bamtofastq -i {params.unmapped_bam} -fq {output.r1_clean} -fq2 {output.r2_clean}
#         """


### Setup Metaphlan. Run Metaphlan on samples to make abundance tables ###

# def metaphlan_merge_inputs(wildcards):
#     files = expand("results/metaphlan_out/{sample}/{sample}.metaphlan_profile.txt",
#         zip, sample=samples["sample"], dataset=samples["dataset"])
#     return files


# rule metaphlan_setup:
#     output:
#         metaphlan_db=directory("resources/metaphlan_db"),
#         metaphlan_db_file="resources/metaphlan_db/{}.rev.1.bt2".format(config["metaphlan_idx"])
#     conda: 
#         "../envs/metaphlan.yml"
#     params:
#         metaphlan_idx = config["metaphlan_idx"] # Index for metaphlan
#     threads: 10
#     resources:
#         mem="50g",
#         time="04:00:00"
#     shell:
#         """
#         metaphlan --install --index {params.metaphlan_idx} --bowtie2db {output.metaphlan_db} --nproc {threads}
#         """


# rule metaphlan:
#     input:
#         metaphlan_db = rules.metaphlan_setup.output.metaphlan_db,
#         r1_clean = "results/bwa_out/{sample}/{sample}.fastp_bwa.r1.fastq",
#         r2_clean = "results/bwa_out/{sample}/{sample}.fastp_bwa.r2.fastq"
#     output:
#         profile = "results/metaphlan_out/{sample}/{sample}.metaphlan_profile.txt",
#         bowtie_out = "results/metaphlan_out/{sample}/{sample}.bowtie2.bz2"
#     conda: 
#         "../envs/metaphlan.yml"
#     params:
#         metaphlan_idx = config["metaphlan_idx"] # Index for metaphlan
#     threads: 20
#     resources:
#         mem="50g",
#         time="04:00:00"
#     shell:
#         """
#         metaphlan {input.r1_clean},{input.r2_clean} \
#         --bowtie2out {output.bowtie_out} \
#         --index {params.metaphlan_idx} \
#         --bowtie2db {input.metaphlan_db} \
#         --nproc {threads} \
#         --input_type fastq \
#         --unknown_estimation \
#         -o {output.profile}
#         """


# rule metaphlan_merge:
#     input:
#         metaphlan_merge_inputs
#     output:
#         "results/metaphlan_merged/merged_metaphlan_profile.tsv"
#     conda:
#         "../envs/metaphlan.yml"
#     shell:
#         """
#         merge_metaphlan_tables.py {input} > {output}
#         """


# rule metaphlan_species_abundance:
#     input:
#         "results/metaphlan_merged/merged_metaphlan_profile.tsv"
#     output:
#         "results/metaphlan_merged/merged_metaphlan_profile_species.tsv"
#     conda:
#         "../envs/metaphlan.yml"
#     shell:
#         """
#         grep -E "s__|clade|UNKNOWN" {input} | sed 's/^.*s__//g' \
#         | cut -f1,3- | sed -e 's/clade_name/sample/g' > {output}
#         """


# rule metaphlan_genus_abundance:
#     input:
#         "results/metaphlan_merged/merged_metaphlan_profile.tsv"
#     output:
#         "results/metaphlan_merged/merged_metaphlan_profile_genus.tsv"
#     conda:
#         "../envs/metaphlan.yml"
#     shell:
#         """
#         grep -E "g__|clade|UNKNOWN" {input} | sed 's/^.*g__//g' \
#         | grep -v s__ |cut -f1,3- | sed -e 's/clade_name/sample/g' > {output}
#         """

# rule metaphlan_unifrac:
#     input:
#         "results/metaphlan_merged/merged_metaphlan_profile.tsv"
#     output:
#         "results/metaphlan_merged/merged_metaphlan_unifrac_matrix.txt"
#     params:
#         "/home/jsj3921/.conda/envs/snakemake/pkgs/metaphlan-3.0.13-pyhb7b1952_0/site-packages/metaphlan/utils/"
#     conda:
#         "../envs/metaphlan.yml"
#     shell:
#         """
#         module load R/4.1.1
#         Rscript {params}calculate_unifrac.R {input} {params}mpa_v30_CHOCOPhlAn_201901_species_tree.nwk {output}
#         """

# rule hclust:
#     input:
#         "results/metaphlan_merged/merged_metaphlan_profile_species.tsv"
#     output:
#         report("results/metaphlan_merged/merged_metaphlan_hclust_species.png", caption="report/hclust.rst", category="METAPHLAN")
#     conda:
#         "../envs/hclust.yml"
#     shell:
#         """
#         hclust2.py -i {input} -o {output} --f_dist_f braycurtis --s_dist_f braycurtis --cell_aspect_ratio 0.5 -l --flabel_size 10 --slabel_size 10 --max_flabel_len 100 --max_slabel_len 100 --minv 0.1 --dpi 300
#         """

# rule hclust_genus:
#     input:
#         "results/metaphlan_merged/merged_metaphlan_profile_genus.tsv"
#     output:
#         report("results/metaphlan_merged/merged_metaphlan_hclust_genus.png", caption="report/hclust_genus.rst", category="METAPHLAN")
#     conda:
#         "../envs/hclust.yml"
#     shell:
#         """
#         hclust2.py -i {input} -o {output} --f_dist_f braycurtis --s_dist_f braycurtis --cell_aspect_ratio 0.5 -l --flabel_size 10 --slabel_size 10 --max_flabel_len 100 --max_slabel_len 100 --minv 0.1 --dpi 300
#         """



# # Single samples assemblies 

# rule megahit:
#     input:
#         r1_clean = "results/bwa_out/{sample}/{sample}.fastp_bwa.r1.fastq",
#         r2_clean = "results/bwa_out/{sample}/{sample}.fastp_bwa.r2.fastq",
#     output:
#         scaffolds = "results/megahit_out/{sample}/final.contigs.fa",
#         outdirec = directory("results/megahit_out/{sample}")
#     threads: 20
#     resources:
#         mem="50g",
#         time="10:00:00"
#     shell:
#         """
#         module load megahit/1.0.6.1
#         megahit -t {threads} -m 0.9 -1 {input.r1_clean} -2 {input.r2_clean} -o {output.outdirec}
#         """

# rule quast:
#     input:
#         scaffolds = "results/megahit_out/{sample}/final.contigs.fa",
#     output:
#         direc=directory("results/quast_out/megahit/{sample}"),
#         report="results/quast_out/megahit/{sample}/report.html"
#     threads: 1
#     conda:
#         "../envs/genome_qc.yml"
#     shell:
#         "quast.py -o {output.direc} --threads {threads} --min-contig 0 -L {input}"

# rule multiqc_quast:
#     input:
#         quast_reports=expand("results/quast_out/megahit/{sample}/report.html", zip, sample=samples["sample"], dataset=samples["dataset"])
#     output:
#         outDir=directory("results/quast_out/megahit/multiqc"),
#         multiqc_report = "results/quast_out/megahit/multiqc/report.html"
#     shell:
#         """
#         module load multiqc
#         multiqc --outdir {output.outDir} {input.quast_reports}
#         """

# rule drop_short_contigs:
#     input:
#         "results/megahit_out/{sample}/final.contigs.fa"
#     output:
#         "results/megahit_out/megahit_g1000/{sample}.megahit_g1000.fa"
#     conda:
#         "../envs/seq_processing.yml"
#     script:
#         "scripts/parse_contigs.py"




# rule concatenate_assemblies:
#     input:
#         expand("results/{dataset}/assembly/megahit_g1000/{sample}.megahit_g1000.fa",
#             zip, sample=samples["sample"], dataset=samples["dataset"])
#     output:
#         "results/allDatasets/single_sample_assemblies/allSamples.megahit_g1000.fa"
#     shell:
#         """
#         cat {input} > {output}
#         """

# rule quast_g1000:
#     input:
#         scaffolds = "results/allDatasets/single_sample_assemblies/allSamples.megahit_g1000.fa"
#     output:
#         direc=directory("results/allDatasets/single_sample_assemblies/quast/monoassemblies_quast"),
#         report="results/allDatasets/single_sample_assemblies/quast/monoassemblies_quast/report.html",
#         tsv_report=report("results/allDatasets/single_sample_assemblies/quast/monoassemblies_quast/report.tsv", caption="report/quast_g1000.rst", category="ASSEMBLY", subcategory="Single Sample"),
#         pdf_report=report("results/allDatasets/single_sample_assemblies/quast/monoassemblies_quast/report.pdf", caption="report/quast_g1000.rst", category="ASSEMBLY", subcategory="Single Sample")

#     threads: 1
#     conda:
#         "../envs/genome_qc.yml"
#     shell:
#         "quast.py -o {output.direc} --threads {threads} --min-contig 0 -L {input}"

        
