# Note that the metaphlan implementation of this does work; however the kaiju
#   implementation is still in alpha-alpha-testing

# Add the following lines to rule all if you want to run this.
        # directory("../resources/kneaddata"), # KneadData DB install
        # expand("../results/{dataset}/kneaddata/{sample}/{sample}_R1_001_kneaddata_paired_1.fastq",
        #     zip, sample=samples["sample"], dataset=samples["dataset"]), # KneadData R2
        # expand("../results/{dataset}/kneaddata/{sample}/{sample}_R1_001_kneaddata_paired_2.fastq",
        #     zip, sample=samples["sample"], dataset=samples["dataset"]), # KneadData R1
        # "../results/allDatasets/metaphlan_kneaddata/merged_abundance_table.KD_allDatasets.txt", #  Kneaddata metaphlan merged table
        # "../results/allDatasets/metaphlan_kneaddata/merged_abundance_table.species.KD_allDatasets.txt",
        # "../results/allDatasets/metaphlan_kneaddata/abundance_heatmap_KD_species.allDatasets.png", # KD Metaphlan hclust
        # "../results/allDatasets/compare/abundance_heatmap_KD_BWA_species.allDatasets.png",
        # expand("../results/{dataset}/abundance/kaiju_refseq_kneaddata/{sample}.kaiju_refseq.txt",
        #     zip, sample=samples["sample"], dataset=samples["dataset"])




### Setup and Execute KneadData ###


rule kneaddata_setup:
    output:
        kneaddata_db = directory("../resources/kneaddata"),
        index = "../resources/kneaddata/hg37dec_v0.1.rev.2.bt2"
    conda:
        "../envs/kneaddata.yml"
    shell:
        """
        kneaddata_database --download human_genome bowtie2 {output.kneaddata_db}
        """


rule kneaddata:
    input:
        kneaddata_db = rules.kneaddata_setup.output.kneaddata_db,
        r1 = get_r1,
        r2 = get_r2
    output:
        outDir = directory("../results/{dataset}/kneaddata/{sample}"),
        cleanR1 = "../results/{dataset}/kneaddata/{sample}/{sample}_R1_001_kneaddata_paired_2.fastq",
        cleanR2 = "../results/{dataset}/kneaddata/{sample}/{sample}_R1_001_kneaddata_paired_1.fastq"
    params:
        db_index = "../resources/kneaddata/hg37dec_v0.1.1.bt2"
    conda:
        "../envs/kneaddata.yml"
    threads: 25
    shell:
        """
        kneaddata \
        --input {input.r1} \
        --input {input.r2} \
        --output {output.outDir} \
        --reference-db {params.db_index} \
        -t {threads} \
        --bypass-trf
        """

rule KD_metaphlan:
    input:
        metaphlan_db = rules.metaphlan_setup.output.metaphlan_db,
        cleanFastQ1 = "../results/{dataset}/kneaddata/{sample}/{sample}_R1_001_kneaddata_paired_1.fastq",
        cleanFastQ2 = "../results/{dataset}/kneaddata/{sample}/{sample}_R1_001_kneaddata_paired_2.fastq"
    output:
        profile = "../results/{dataset}/abundance/metaphlan_kneaddata/{sample}.metaphlan_profile_KD.txt",
        bowtie_out = "../results/{dataset}/abundance/metaphlan_kneaddata/{sample}.bowtie2.bz2"
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
        --unknown_estimation \
        -o {output.profile}
        """

def metaphlan_merge_kneaddata_inputs(wildcards):
    files = expand("../results/{dataset}/abundance/metaphlan_kneaddata/{sample}.metaphlan_profile_KD.txt",
        zip, sample=samples["sample"], dataset=samples["dataset"])
    return files


rule KD_metaphlan_merge:
    input:
        metaphlan_merge_kneaddata_inputs
    output:
        "../results/allDatasets/metaphlan_kneaddata/merged_abundance_table.KD_allDatasets.txt"
    conda:
        "../envs/metaphlan.yml"
    shell:
        """
        merge_metaphlan_tables.py {input} > {output}
        """


rule KD_metaphlan_species_abundance:
    input:
        "../results/allDatasets/metaphlan_kneaddata/merged_abundance_table.KD_allDatasets.txt"
    output:
        "../results/allDatasets/metaphlan_kneaddata/merged_abundance_table.species.KD_allDatasets.txt"
    conda:
        "../envs/metaphlan.yml"
    shell:
        """
        grep -E "s__|clade|UNKNOWN" {input} | sed 's/^.*s__//g' \
        | cut -f1,3- | sed -e 's/clade_name/sample/g' > {output}
        """


rule KD_hclust:
    input:
        "../results/allDatasets/metaphlan_kneaddata/merged_abundance_table.species.KD_allDatasets.txt"
    output:
        "../results/allDatasets/metaphlan_kneaddata/abundance_heatmap_KD_species.allDatasets.png"
    conda:
        "../envs/hclust.yml"
    shell:
        """
        hclust2.py -i {input} -o {output} --f_dist_f braycurtis --s_dist_f braycurtis --cell_aspect_ratio 0.5 -l --flabel_size 6 --slabel_size 15 --max_flabel_len 100 --max_slabel_len 100 --minv 0.1 --dpi 300
        """

rule KDBWA_metaphlan_merge:
    input:
        KD = metaphlan_merge_kneaddata_inputs,
        BWA = metaphlan_merge_inputs
    output:
        "../results/allDatasets/compare/merged_abundance_table.KD_BWA_allDatasets.txt"
    conda:
        "../envs/metaphlan.yml"
    shell:
        """
        merge_metaphlan_tables.py {input.KD} {input.BWA} > {output}
        """

rule KDBWA_metaphlan_species_abundance:
    input:
        "../results/allDatasets/compare/merged_abundance_table.KD_BWA_allDatasets.txt"
    output:
        "../results/allDatasets/compare/merged_abundance_table.species.KD_BWA_allDatasets.txt"
    conda:
        "../envs/metaphlan.yml"
    shell:
        """
        grep -E "s__|clade|UNKNOWN" {input} | sed 's/^.*s__//g' \
        | cut -f1,3- | sed -e 's/clade_name/sample/g' > {output}
        """



rule KDBWA_hclust:
    input:
        "../results/allDatasets/compare/merged_abundance_table.species.KD_BWA_allDatasets.txt"
    output:
        "../results/allDatasets/compare/abundance_heatmap_KD_BWA_species.allDatasets.png"
    conda:
        "../envs/hclust.yml"
    shell:
        """
        hclust2.py -i {input} -o {output} --f_dist_f braycurtis --s_dist_f braycurtis --cell_aspect_ratio 1.0 -l --flabel_size 8 --slabel_size 10 --max_flabel_len 100 --max_slabel_len 100 --minv 0.1 --dpi 300
        """
##

rule KD_sort:
    input:
        cleanFastQ1 = "../results/{dataset}/kneaddata/{sample}/{sample}_R1_001_kneaddata_paired_1.fastq",
        cleanFastQ2 = "../results/{dataset}/kneaddata/{sample}/{sample}_R1_001_kneaddata_paired_2.fastq"
    output:
        sortedR1 = "../results/{dataset}/kneaddata_sorted/{sample}_R1_001_kneaddata_paired_1.sorted.fastq",
        sortedR2 = "../results/{dataset}/kneaddata_sorted/{sample}_R1_001_kneaddata_paired_2.sorted.fastq"
    shell:
        """
        cat {input.cleanFastQ1} | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > {output.sortedR1}
        cat {input.cleanFastQ2} | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > {output.sortedR2}
        """


rule kaiju_refseq_KD:
    input:
        kaiju_sb = rules.kaiju_db.output.tar,
        cleanFastQ1 = "../results/{dataset}/kneaddata/{sample}/{sample}_R1_001_kneaddata_paired_1.fastq",
        cleanFastQ2 = "../results/{dataset}/kneaddata/{sample}/{sample}_R1_001_kneaddata_paired_2.fastq"
    output:
        profile = "../results/{dataset}/abundance/kaiju_refseq_kneaddata/{sample}.kaiju_refseq.txt",
    params:
        db_path = "../resources/kaiju_head/kaijuDB",
        mode = "mem"
    threads: 25
    shell:
        """
        ../resources/kaiju_head/kaijuDir/kaiju -z {threads} \
        -t {params.db_path}/nodes.dmp \
        -f {params.db_path}/kaiju_db_refseq.fmi \
        -i {input.cleanFastQ1} \
        -j {input.cleanFastQ2} \
        -a {params.mode} \
        -o {output.profile}
        """

# rule KD_kaiju_refseq:
#     input:
#         kaiju_sb = rules.kaiju_db.output.tar,
#         sortedR1 = "../results/{dataset}/kneaddata_sorted/{sample}_R1_001_kneaddata_paired_1.sorted.fastq",
#         sortedR2 = "../results/{dataset}/kneaddata_sorted/{sample}_R1_001_kneaddata_paired_2.sorted.fastq"
#     output:
#         profile = "../results/{dataset}/abundance/kaiju_refseq_KD/{sample}.kaiju_refseq_KD.txt"
#     params:
#         db_path = "../resources/kaiju_head/kaijuDB",
#         mode = "mem"
#     threads: 25
#     shell:
#         """
#         ../resources/kaiju_head/kaijuDir/kaiju -z {threads} \
#         -t {params.db_path}/nodes.dmp \
#         -f {params.db_path}/kaiju_db_refseq.fmi \
#         -i {input.sortedR1} \
#         -j {input.sortedR2} \
#         -a {params.mode} \
#         -o {output.profile}
#         """


def kaiju_merge_inputs_KD(wildcards):
    files = expand("../results/{dataset}/abundance/kaiju_refseq_KD/{sample}.kaiju_refseq_KD.txt",
        zip, sample=samples["sample"], dataset=samples["dataset"])
    return files


rule KD_kaiju_merge_KD:
    input: 
        kaiju_merge_inputs_KD
    output:
        "../results/allDatasets/kaiju_KD/kaiju_abundance_table_KD.species.allDatasets.txt"
    params:
        db_path = "../resources/kaiju_head/kaijuDB",
        taxa = "species"
    shell:
        """
        ../resources/kaiju_head/kaijuDir/kaiju2table \
        -t {params.db_path}/nodes.dmp \
        -n {params.db_path}/names.dmp \
        -r {params.taxa} \
        -o {output} \
        {input}
        """
 
