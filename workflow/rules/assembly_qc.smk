
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

