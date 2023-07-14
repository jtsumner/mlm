
############################
###   PART 1A: MEGAHIT   ###
############################

rule megahit:
    input:
        r1_clean = get_final_read1,
        r2_clean = get_final_read2
    output:
        scaffolds = "results/megahit_out/{sample}/{sample}.contigs.fa"
    params:
        out_dir = "results/megahit_out/{sample}"
    threads: 20
    resources:
        mem="100g",
        time="10:00:00"
    shell:
        """
        module load megahit/1.0.6.1
        megahit -t {threads} -m 0.9 -1 {input.r1_clean} -2 {input.r2_clean} --out-prefix {wildcards.sample} -o {params.out_dir}_tmp
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

    threads: 25
    resources:
        mem="100g",
        time="10:00:00"
    shell:
        """
        module load spades/3.14.1
        spades.py -1 {input.r1_clean} -2 {input.r2_clean} -o {params.out_dir} -t {threads} -m 100 --meta
        """
# -k 21,33,55,77,99,127 --only-assembler

############################
###   PART 2: FILTER     ###
############################

rule drop_short_contigs_megahit:
    input:
        "results/megahit_out/{sample}/{sample}.contigs.fa"
    output:
        "results/megahit_parsed/{sample}/{sample}.parsed_assembly.fa"
    conda:
        "../envs/biopython.yml"
    script:
        "../scripts/parse_contigs.py"

rule drop_short_contigs_spades:
    input:
        scaffolds = "results/spades_out/{sample}/scaffolds.fasta"
    output:
        parsed = "results/spades_parsed/{sample}/{sample}.fa"
    params:
        min_length = "1000",
        assembler="spades"
    conda:
        "../envs/biopython.yml"
    shell:
        """
        python3 workflow/scripts/parse_contigs.py --sample {wildcards.sample} \
            --scaffolds {input.scaffolds} \
            --parsed {output.parsed} \
            --min_length {params.min_length} \
            --assembler {params.assembler}
        """

############################
###    PART 3: QUAST     ###
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
        multiqc --outdir {params.out_dir} --dirs --dirs-depth 2 results/quast_out/
        """
