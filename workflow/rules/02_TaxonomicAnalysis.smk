
############################
### PART 1A: METAPHLAN3  ###
############################

### Setup Metaphlan. Run Metaphlan on samples to make abundance tables ###
#mpa_vOct22_CHOCOPhlAnSGB_202212
def metaphlan_merge_inputs(wildcards):
    files = expand("results/metaphlan_out/{sample}/{sample}.metaphlan_profile.txt",
        sample=samples["sample"])
    return files

#        metaphlan_db_file="resources/metaphlan_db/{}.rev.1.bt2".format(config["metaphlan_idx"])

rule metaphlan_setup:
    output:
        metaphlan_db=directory("resources/metaphlan_db"),
        metaphlan_db_file="resources/metaphlan_db/{}.rev.2.bt2l".format(config["metaphlan_idx"])

    conda: 
        "../envs/metaphlan.yml"
    params:
        metaphlan_idx = config["metaphlan_idx"] # Index for metaphlan
    threads: 10
    resources:
        mem="10g",
        time="04:00:00"
    shell:
        """
        metaphlan --install --index {params.metaphlan_idx} --bowtie2db {output.metaphlan_db} --nproc {threads}
        """

rule metaphlan:
    input:
        metaphlan_db = rules.metaphlan_setup.output.metaphlan_db,
        r1_clean = "results/bwa_out/{sample}/{sample}.fastp_bwa.r1.fastq",
        r2_clean = "results/bwa_out/{sample}/{sample}.fastp_bwa.r2.fastq"
    output:
        profile = "results/metaphlan_out/{sample}/{sample}.metaphlan_profile.txt",
        bowtie_out = "results/metaphlan_out/{sample}/{sample}.bowtie2.bz2"
    conda: 
        "../envs/metaphlan.yml"
    params:
        metaphlan_idx = config["metaphlan_idx"] # Index for metaphlan
    threads: 10
    resources:
        mem="20g",
        time="04:00:00"
    shell:
        """
        metaphlan {input.r1_clean},{input.r2_clean} \
        --bowtie2out {output.bowtie_out} \
        --index {params.metaphlan_idx} \
        --bowtie2db {input.metaphlan_db} \
        --nproc {threads} \
        --input_type fastq \
        --unclassified_estimation \
        -t rel_ab_w_read_stats \
        -o {output.profile}
        """

rule metaphlan_merge:
    input:
        metaphlan_merge_inputs
    output:
        "results/metaphlan_merged/merged_metaphlan_profile.tsv"
    shell:
        """
        python3 workflow/scripts/merge_metaphlan_tables_abs.py {input} > {output}
        """

rule metaphlan_species_abundance:
    input:
        "results/metaphlan_merged/merged_metaphlan_profile.tsv"
    output:
        "results/metaphlan_merged/merged_metaphlan_profile_species.tsv"
    conda:
        "../envs/metaphlan.yml"
    shell:
        """
        grep -E "s__|clade|UNKNOWN" {input} | sed 's/^.*s__//g' \
        | cut -f1,3- | sed -e 's/clade_name/sample/g' > {output}
        """

rule metaphlan_genus_abundance:
    input:
        "results/metaphlan_merged/merged_metaphlan_profile.tsv"
    output:
        "results/metaphlan_merged/merged_metaphlan_profile_genus.tsv"
    conda:
        "../envs/metaphlan.yml"
    shell:
        """
        grep -E "g__|clade|UNKNOWN" {input} | sed 's/^.*g__//g' \
        | grep -v s__ |cut -f1,3- | sed -e 's/clade_name/sample/g' > {output}
        """

rule metaphlan_unifrac:
    input:
        "results/metaphlan_merged/merged_metaphlan_profile.tsv"
    output:
        "results/metaphlan_merged/merged_metaphlan_unifrac_matrix.txt"
    params:
        "/home/jsj3921/.conda/envs/snakemake/pkgs/metaphlan-3.0.13-pyhb7b1952_0/site-packages/metaphlan/utils/"
    conda:
        "../envs/metaphlan.yml"
    shell:
        """
        module load R/4.1.1
        Rscript {params}calculate_unifrac.R {input} {params}mpa_v30_CHOCOPhlAn_201901_species_tree.nwk {output}
        """

rule hclust:
    input:
        "results/metaphlan_merged/merged_metaphlan_profile_species.tsv"
    output:
        report("results/metaphlan_merged/merged_metaphlan_hclust_species.png", caption="report/hclust.rst", category="METAPHLAN")
    conda:
        "../envs/hclust.yml"
    shell:
        """
        hclust2.py -i {input} -o {output} --f_dist_f braycurtis --s_dist_f braycurtis --cell_aspect_ratio 0.5 -l --flabel_size 10 --slabel_size 10 --max_flabel_len 100 --max_slabel_len 100 --minv 0.1 --dpi 300
        """

rule hclust_genus:
    input:
        "results/metaphlan_merged/merged_metaphlan_profile_genus.tsv"
    output:
        report("results/metaphlan_merged/merged_metaphlan_hclust_genus.png", caption="report/hclust_genus.rst", category="METAPHLAN")
    conda:
        "../envs/hclust.yml"
    shell:
        """
        hclust2.py -i {input} -o {output} --f_dist_f braycurtis --s_dist_f braycurtis --cell_aspect_ratio 0.5 -l --flabel_size 10 --slabel_size 10 --max_flabel_len 100 --max_slabel_len 100 --minv 0.1 --dpi 300
        """

use rule metaphlan as metaphlan_bowtie with:
    input:
        metaphlan_db = rules.metaphlan_setup.output.metaphlan_db,
        r1_clean = get_final_read1,
        r2_clean = get_final_read2
    output:
        profile = "results/metaphlan_bowtie_out/{sample}/{sample}.metaphlan_profile.txt",
        bowtie_out = "results/metaphlan_bowtie_out/{sample}/{sample}.bowtie2.bz2"

use rule metaphlan_merge as metaphlan_merge_bowtie with:
    input:
        expand("results/metaphlan_bowtie_out/{sample}/{sample}.metaphlan_profile.txt", zip, sample=samples["sample"], dataset=samples["dataset"])
    output:
        "results/metaphlan_bowtie_out/merged_metaphlan_profile.tsv"

use rule metaphlan_genus_abundance as metaphlan_bowtie_genus_abundance with:
    input:
        "results/metaphlan_bowtie_out/merged_metaphlan_profile.tsv"
    output:
        "results/metaphlan_bowtie_out/merged_metaphlan_profile_genus.tsv"

use rule metaphlan_species_abundance as metaphlan_bowtie_species_abundance with:
    input:
        "results/metaphlan_bowtie_out/merged_metaphlan_profile.tsv"
    output:
        "results/metaphlan_bowtie_out/merged_metaphlan_profile_species.tsv"

rule metaphlan_bbmerge:
    input:
        metaphlan_db = rules.metaphlan_setup.output.metaphlan_db,
        merged_reads = get_final_merged_read #get_final_read1
    output:
        profile = "results/metaphlan_bbmerge_out/{sample}/{sample}.metaphlan_profile.txt",
        bowtie_out = "results/metaphlan_bbmerge_out/{sample}/{sample}.bowtie2.bz2"
    conda: 
        "../envs/metaphlan.yml"
    params:
        metaphlan_idx = config["metaphlan_idx"] # Index for metaphlan
    threads: 10
    resources:
        mem="20g",
        time="04:00:00"
    shell:
        """
        metaphlan {input.merged_reads} \
        --bowtie2out {output.bowtie_out} \
        --index {params.metaphlan_idx} \
        --bowtie2db {input.metaphlan_db} \
        --nproc {threads} \
        --input_type fastq \
        --unclassified_estimation \
        -t rel_ab_w_read_stats \
        -o {output.profile}
        """
############################
###  PART 1B: KRACKEN2   ###
############################

rule kraken2:
    """
    Performs taxnomic classification with Kraken2 

    Outputs a kraken2-style report and metaphlan-style report with a
    script from KrakenTools
    """
    input: 
        r1_clean = get_final_read1,
        r2_clean = get_final_read2
    output:
        kraken_out = "results/kraken/{sample}/{sample}_kraken2out.txt",
        kraken_report = "results/kraken/{sample}/{sample}_kraken2report.txt"
    threads: 20
    resources:
        mem="650G",
        time="01:00:00",
        partition="genomics-himem"
    params:
        kraken_db = "/projects/b1188/bmo/mlm/resources/kraken_db" #/software/kraken/database/kraken_db
    shell:
        """
        module load kraken/2
        kraken2 --threads {threads} \
            --output {output.kraken_out} \
            --report {output.kraken_report} \
            --confidence 0.7 \
            --gzip-compressed \
            --minimum-hit-groups 3 \
            --db {params.kraken_db} \
            --paired {input.r1_clean} {input.r2_clean}
        """ 

rule kraken_mpa:
    """
    Outputs a metaphlan-style report with a script from KrakenTools
    """
    input: 
        kraken_report = "results/kraken/{sample}/{sample}_kraken2report.txt"
    output:
        kraken_mpa = "results/kraken/{sample}/{sample}_mpa.tsv"
    threads: 1
    resources:
        mem="3G",
        time="0:10:00"
    shell:
        """
        workflow/scripts/kreport2mpa.py -r {input.kraken_report} \
            -o {output.kraken_mpa} \
            --display-header \
            --no-intermediate-ranks
        """ 

rule merge_kraken:
    """
    Outputs a metaphlan-style report with a script from KrakenTools
    """
    input: 
        kraken_reports = expand("results/kraken/{sample}/{sample}_kraken2report.txt", sample=samples["sample"]),
        kraken_mpas = expand("results/kraken/{sample}/{sample}_mpa.tsv", sample=samples["sample"])
    output:
        merged_reports = "results/kraken/merged_kraken_report_profile.tsv",
        merged_mpa = "results/kraken/merged_kraken_mpa_profile.tsv"
    threads: 1
    resources:
        mem="15G",
        time="0:10:00"
    shell:
        """
        workflow/scripts/combine_kreports.py -r {input.kraken_reports} \
            -o {output.merged_reports} \
            --display-headers
        
        workflow/scripts/combine_mpa.py -i {input.kraken_mpas} \
            -o {output.merged_mpa} \
        """ 

rule bracken:
    """
    Performs abundance estimation from with Kraken2 classification
    """
    input: 
        kraken_out = "results/kraken/{sample}/{sample}_kraken2out.txt",
        kraken_report = "results/kraken/{sample}/{sample}_kraken2report.txt"
    output:
        bracken_out = "bracken_out/{sample}/{sample}.bracken",
        bracken_report = "bracken_out/{sample}/{sample}.breport"
    threads: 2
    conda:
        "../envs/bracken.yml"
    resources:
        mem="5G",
        time="0:10:00"
    params:
        read_length = "100",
        taxonomic_level = "S",
        read_threshold = "10",
        kraken_db = "/projects/b1188/bmo/mlm/resources/kraken_db" #/software/kraken/database/kraken_db
    shell:
        """
        module load kraken/2
        bracken-build -d krakendb -t 8 -k 35 -l 100
        bracken -d {params.kraken_db} \
            -i {input.kraken_report} \
            -r {params.read_length} \
            -l {params.taxonomic_level} \
            -t {params.read_threshold} \
            -o {output.bracken_out} \
            -w {output.bracken_report} \
        """  

############################
###   PART 1C: METAXA2   ###
############################

rule metaxa2:
    input:
        r1_clean = get_final_read1,
        r2_clean = get_final_read2
    output:
        out = "results/metaxa2/{sample}/{sample}_metaxa2.taxonomy.txt"
    params:
        base_out = "results/metaxa2/{sample}/{sample}_metaxa2"
    threads: 25
    resources:
        mem="30G",
    shell:
        """
        module load hmmer/3.1b2 blast/2.7.1 mafft/7.407 python/anaconda3 metaxa2/2.2

        metaxa2 -1 {input.r1_clean} \
            -2 {input.r1_clean} \
            --mode metagenome \
            -f fastq \
            -z gzip \
            -g ssu \
            -p /software/metaxa2/2.2/metaxa2_db/SSU/HMMs/ \
            -o {params.base_out} \
            --cpu 24 \
            --multi_thread T \
            --plus T \
            --graphical F \
            --fasta F \

        """

#  metaxa2_ttt -i 20221221-Comunal-Zymo_metaxa2.taxonomy.txt -o test -t A,b 
# --cpu 24 --multi_thread T --unknown T -r 1 --distance=0


############################
###    PART 2: HUMANN    ###
############################

rule humann:
    input:
        metaphlan_db = rules.metaphlan_setup.output.metaphlan_db,
        merged_reads = get_final_merged_read #get_final_read1
    output:
        gene_fam = "results/humann_out/{sample}/{sample}_genefamilies.tsv",
        path_cov = "results/humann_out/{sample}/{sample}_pathcoverage.tsv",
        path_abund = "results/humann_out/{sample}/{sample}_pathabundance.tsv"
    conda: 
        "../envs/humann.yml"
    params:
        metaphlan_idx = config["metaphlan_idx"] # Index for metaphlan
    threads: 20
    resources:
        mem="60G",
        time="02:00:00"
    shell:
        """
        outdir=$(dirname {output.gene_fam})
        humann --input {input.merged_reads} \
            --output $outdir \
            --threads {threads} \
            --nucleotide-database /projects/b1180/software/conda_envs/humann/lib/python3.7/site-packages/humann/data/chocophlan/ \
            --protein-database /projects/b1180/software/conda_envs/humann/lib/python3.7/site-packages/humann/data/uniref/ \
            --metaphlan-options="--index {params.metaphlan_idx} --bowtie2db {input.metaphlan_db}"
        """