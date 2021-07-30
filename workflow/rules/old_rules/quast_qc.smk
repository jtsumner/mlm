
rule quast_raw:
    input:
        "../data/spades_assemblies/{sample}/scaffolds.fasta"
    output:
        direc=directory("../data/quast/raw_scaffolds/results_{sample}"),
        report="../data/quast/raw_scaffolds/results_{sample}/report.html"
    threads: 1
    conda:
        "../envs/genome_qc.yml"
    shell:
        "quast.py -o {output.direc} --threads {threads} {input}"


rule quast_parsed:
    input:
        "../data/parsed_assemblies/{sample}_greater1kb_scaffolds.fasta"
    output:
        direc=directory("../data/quast/parsed_scaffolds/results_{sample}"),
        report="../data/quast/parsed_scaffolds/results_{sample}/report.html"
    threads: 1
    conda:
        "../envs/genome_qc.yml"
    shell:
        "quast.py -o {output.direc} --threads {threads} {input}"