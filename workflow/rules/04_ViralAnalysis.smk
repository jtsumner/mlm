
#### Download Vibrant database for the prediction of viral contigs ####
rule vibrant_download_db:
    output:
        hack="resources/vibrant/vibrant_db_downloaded.txt"
    threads: 2
    resources:
        mem="25g"
    conda:
        "../envs/vibrant.yml"
    shell:
        """
        download-db.sh
        touch resources/vibrant/vibrant_db_downloaded.txt
        """

def vibrant_input(wildcards):
    files = expand("results/spades_parsed/{sample}/{sample}.fa", sample=samples["sample"])
    return files

#### Prepare Spades output for VContact2 ####

rule concat_parsed_assemblies:
    input:
        files = vibrant_input
    output:
        temp("results/vibrant_out/parsed_scaffolds.fa")
    shell:
        "cat {input.files} > {output}"


#### Run Vibrant on parsed assemblies (i.e., contigs > 1kb ) to predict viral contigs ####
rule vibrant_parsed:
    input:
        fasta="results/vibrant_out/parsed_scaffolds.fa",
        hack="resources/vibrant/vibrant_db_downloaded.txt"
    output:
        #direc=directory("results/vibrant_output/{sample}/VIBRANT_{sample}"),
        #phages="results/vibrant_output/{sample}/VIBRANT_{sample}/VIBRANT_phages_{sample}/{sample}.phages_combined.faa"
        direc=directory("results/vibrant_output/VIBRANT_parsed_scaffolds"),
        phages="results/vibrant_output/VIBRANT_parsed_scaffolds/VIBRANT_phages_parsed_scaffolds/parsed_scaffolds.phages_combined.faa"
    threads: 60
    resources:
        mem="200g",
        time="10:00:00"
    conda:
        "../envs/vibrant.yml"
    shell:
         "VIBRANT_run.py -i {input.fasta} -d $VIBRANT_DATA_PATH/databases -folder results/vibrant_output/VIBRANT_parsed_scaffolds/ -t {threads}"

#### Prepare Vibrant output for VContact2 ####



rule simplify_vibrant_AA:
    input:
        "results/vibrant_output/VIBRANT_parsed_scaffolds/VIBRANT_phages_parsed_scaffolds/parsed_scaffolds.phages_combined.faa"
    output:
        "results/vibrant_output/VIBRANT_parsed_scaffolds/VIBRANT_phages_parsed_scaffolds/parsed_scaffolds.phages_combined.simple.faa"
    conda:
        "../envs/vibrant.yml"
    script:
        "../scripts/simplify_faa-ffn.py"


rule generate_gene2genome:
    input:
        "results/vibrant_output/VIBRANT_parsed_scaffolds/VIBRANT_phages_parsed_scaffolds/parsed_scaffolds.phages_combined.simple.faa"
    output:
        "results/vcontact2_data/g2g_concat_phages_combined.simple.csv"
    conda:
        "../envs/vcontact2.yml"
    shell:
        "vcontact2_gene2genome -p {input} -o {output} -s 'Prodigal-FAA'"

        #wget http://www.paccanarolab.org/static_content/clusterone/cluster_one-1.0.jar

# Pandas 1.0.5, scipy 1.6.0
rule vcontact_parsed:
    input:
        g2gome = "results/vcontact2_data/g2g_concat_phages_combined.simple.csv",
        faa =  "results/vibrant_output/VIBRANT_parsed_scaffolds/VIBRANT_phages_parsed_scaffolds/parsed_scaffolds.phages_combined.simple.faa"
    output:
        v_direc = directory("results/vcontact2_data/vcontact2_output"),
        genome_by_genome = "results/vcontact2_data/vcontact2_output/genome_by_genome_overview.csv"
    threads: 60
    resources:
        mem="80g",
        time="10:00:00",
        nodes=1
    conda:
        "../envs/vcontact2.yml"
    shell:
        """
        vcontact2 --raw-proteins {input.faa} \
            --rel-mode 'Diamond' \
            --proteins-fp {input.g2gome} \
            --db 'ProkaryoticViralRefSeq211-Merged' \
            --pcs-mode MCL \
            --vcs-mode ClusterONE \
            --c1-bin resources/cluster_one/cluster_one-1.0.jar \
            --output-dir results/vcontact2_data/vcontact2_output \
            -t {threads}
        """

