
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
        -t rel_ab \
        -o {output.profile}
        """

rule metaphlan_merge:
    input:
        metaphlan_merge_inputs
    output:
        "results/metaphlan_merged/merged_metaphlan_profile.tsv"
    conda: 
        "../envs/metaphlan.yml"
    shell:
        """
        merge_metaphlan_tables.py {input} > {output}
        """
#       python3 workflow/scripts/merge_metaphlan_tables_abs.py {input} > {output} # for stats based

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
        -t rel_ab \
        -o {output.profile}
        """
        #-t rel_ab_w_read_stats \
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
        time="0:02:00"
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
    params:
        sample_list = list(samples["sample"])
    threads: 1
    resources:
        mem="15G",
        time="0:10:00"
    shell:
        """
        workflow/scripts/combine_kreports.py -r {input.kraken_reports} \
            -o {output.merged_reports} \
            --sample_list {params.sample_list} \
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
        bracken_out = "results/bracken_out/{sample}/{sample}.bracken",
        bracken_report = "results/bracken_out/{sample}/{sample}.breport"
    threads: 1
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
        bracken -d {params.kraken_db} \
            -i {input.kraken_report} \
            -r {params.read_length} \
            -l {params.taxonomic_level} \
            -t {params.read_threshold} \
            -o {output.bracken_out} \
            -w {output.bracken_report} \
        """  

use rule kraken_mpa as bracken_mpa with:
    input: 
        kraken_report = "results/bracken_out/{sample}/{sample}.breport"
    output:
        kraken_mpa = "results/bracken_out/{sample}/{sample}_mpa.tsv"

use rule merge_kraken as merge_bracken with:
    input: 
        kraken_reports = expand("results/bracken_out/{sample}/{sample}.breport", sample=samples["sample"]),
        kraken_mpas = expand("results/bracken_out/{sample}/{sample}_mpa.tsv", sample=samples["sample"])
    output:
        merged_reports = "results/bracken_out/merged_bracken_report_profile.tsv",
        merged_mpa = "results/bracken_out/merged_bracken_mpa_profile.tsv"
#"results/bracken_out/{sample}/{sample}.bracken"
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
        r1_clean = get_final_read1,
        r2_clean = get_final_read2
    output:
        concatenated_reads = temp("results/concat_reads_humann/{sample}.fastq.gz"),
        gene_fam = "results/humann_out/{sample}/{sample}_genefamilies.tsv",
        path_cov = "results/humann_out/{sample}/{sample}_pathcoverage.tsv",
        path_abund = "results/humann_out/{sample}/{sample}_pathabundance.tsv"
    conda: 
        "../envs/humann.yml"
    params:
        metaphlan_idx = config["metaphlan_idx"] # Index for metaphlan
    threads: 20
    resources:
        mem="30G",
        time="03:30:00"
    shell:
        """
        cat {input.r1_clean} {input.r2_clean} > {output.concatenated_reads}

        outdir=$(dirname {output.gene_fam})
        
        humann --input {output.concatenated_reads} \
            --output $outdir \
            --threads {threads} \
            --translated-subject-coverage-threshold 0.0 --nucleotide-subject-coverage-threshold 0.0 --bowtie-options="--very-sensitive-local" \
            --nucleotide-database /projects/b1180/software/conda_envs/humann/lib/python3.7/site-packages/humann/data/chocophlan/ \
            --protein-database /projects/b1180/software/conda_envs/humann/lib/python3.7/site-packages/humann/data/uniref/ \
            --metaphlan-options="-t rel_ab --index {params.metaphlan_idx} --bowtie2db {input.metaphlan_db}"
        """

rule merge_humann:
    input:
        gene_fam = expand("results/humann_out/{sample}/{sample}_genefamilies.tsv", sample=samples["sample"]),
        path_cov = expand("results/humann_out/{sample}/{sample}_pathcoverage.tsv", sample=samples["sample"]),
        path_abund = expand("results/humann_out/{sample}/{sample}_pathabundance.tsv", sample=samples["sample"])
    output:
        gene_fam = "results/humann_out/merged_genefamilies.tsv",
        path_cov = "results/humann_out/merged_pathcoverage.tsv",
        path_abund = "results/humann_out/merged_pathabundance.tsv"
    params:
        humann_dir = "results/humann_out"
    resources:
        time="0:15:00",
        mem="20G",
    threads: 1
    conda: 
        "../envs/humann.yml"
    shell:
        """
        humann_join_tables -i {params.humann_dir} \
            -o {output.gene_fam} \
            --file_name genefamilies \
            --search-subdirectories 

        humann_join_tables -i {params.humann_dir} \
            -o {output.path_cov} \
            --file_name pathcoverage \
            --search-subdirectories

        humann_join_tables -i {params.humann_dir} \
            -o {output.path_abund} \
            --file_name pathabundance \
            --search-subdirectories
        """

rule renorm_humann:
    input:
        gene_fam = "results/humann_out/merged_genefamilies.tsv",
        path_cov = "results/humann_out/merged_pathcoverage.tsv",
        path_abund = "results/humann_out/merged_pathabundance.tsv"
    output:
        gene_fam = "results/humann_out/merged_genefamilies-cpm.tsv",
    params:
        humann_dir = "results/humann_out"
    resources:
        time="01:00:00",
        mem = "200G"
    threads: 1
    conda: 
        "../envs/humann.yml"
    shell:
        """
        humann_renorm_table --input {input.gene_fam} \
            --output {output.gene_fam} \
            --units cpm \
            --update-snames
        """

use rule renorm_humann as renorm_humann_path with:
    input:
        gene_fam = "results/humann_out/merged_pathabundance.tsv",
        path_cov = "results/humann_out/merged_pathcoverage.tsv",
        path_abund = "results/humann_out/merged_pathabundance.tsv"
    output:
        gene_fam = "results/humann_out/merged_pathabundance-cpm.tsv",
    resources:
        time="0:30:00",
        mem = "20G"


rule regroup_humann:
    input:
        gene_fam = "results/humann_out/merged_genefamilies-cpm.tsv",
    output:
        ko_table = "results/humann_out/ko_genefamilies-cpm.tsv",
    params:
        mapping_file = "/projects/b1180/software/conda_envs/humann/lib/python3.7/site-packages/humann/data/utility_mapping/map_ko_uniref90.txt.gz"
    resources:
        time="01:00:00",
        mem = "150G"
    threads: 1
    conda: 
        "../envs/humann.yml",
    shell:
        """
        humann_regroup_table --input {input.gene_fam} \
            -c {params.mapping_file} \
            --output {output.ko_table}
        """

    # rule humann_stratified

    # rule humann_rename_keggs


rule stratify_humann_kegg:
    input:
        kegg = "results/humann_out/ko_genefamilies-cpm.tsv",
    output:
        strat = "results/humann_out/ko_genefamilies-cpm_stratified.tsv",
        unstrat = "results/humann_out/ko_genefamilies-cpm_unstratified.tsv"
    params:
        humann_dir = "results/humann_out"
    resources:
        time="01:00:00",
        mem = "30G"
    threads: 1
    conda: 
        "../envs/humann.yml"
    shell:
        """
        humann_split_stratified_table --input {input.kegg} \
            --output {params.humann_dir} 
        """
#humann_rename_table -i results/humann_out/ko_genefamilies-cpm_unstratified.tsv -n kegg-orthology -o results/humann_out/ko_genefamilies-cpm_unstratified_named.tsv 