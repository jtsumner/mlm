
rule quast_megahit:
    input:
        scaffolds = "results/megahit_out/{sample}/{sample}.contigs.fa"
    output:
        out_dir=directory("results/quast_out/megahit/{sample}"),
        report="results/quast_out/megahit/{sample}/report.html"
    threads: 1
    conda:
        "../envs/genome_qc.yml"
    shell:
        "quast.py -o {output.out_dir} --threads {threads} --min-contig 0 -L {input}"


rule quast_spades:
    input:
        scaffolds="results/spades_out/{sample}/scaffolds.fasta"
    output:
        out_dir=directory("results/quast_out/spades/{sample}"),
        report="results/quast_out/quast/{sample}/report.html"
    threads: 1
    conda:
        "../envs/genome_qc.yml"
    shell:
        "quast.py -o {output.out_dir} --threads {threads} --min-contig 0 -L {input}"


rule multiqc_quast:
    input:
        quast_reports=expand(
            "results/quast_out/megahit/{sample}/report.html", 
            sample=samples["sample"]
        )
    output:
        out_dir=directory("results/quast_out/megahit/multiqc"),
        multiqc_report = "results/quast_out/megahit/multiqc/multiqc_report.html"
    shell:
        """
        module load multiqc
        multiqc --outdir {output.out_dir} results/quast_out/
        """

rule drop_short_contigs:
    input:
        "results/megahit_out/{sample}/{sample}.contigs.fa"
    output:
        "results/megahit_out/megahit_g1000/{sample}.megahit_g1000.fa"
    conda:
        "../envs/seq_processing.yml"
    script:
        "scripts/parse_contigs.py"
