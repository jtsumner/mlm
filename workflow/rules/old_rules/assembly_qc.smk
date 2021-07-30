rule parse_assembly:
    input:
        "../data/spades_assemblies/{sample}/scaffolds.fasta"
    output:
        "../data/parsed_assemblies/{sample}_greater1kb_scaffolds.fasta"
    conda:
        "../envs/seq_processing.yml"
    script:
        "../scripts/parse_contigs.py"

