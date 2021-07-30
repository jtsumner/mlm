
#### Rule to assemble concatenated sample reads with Spades ####
rule assemble:
    input:
        r1="../data/concat_trimmed/{sample}.r1.fastq.gz",
        r2="../data/concat_trimmed/{sample}.r2.fastq.gz"
    output:
        direc=directory("../data/spades_assemblies/{sample}"),
        scaffolds="../data/spades_assemblies/{sample}/scaffolds.fasta"
    threads: 40
    shell:
        "module load spades; spades.py -1 {input.r1} -2 {input.r2} -o {output.direc} -t {threads} -m 100 --meta"

#TODO fix unit so that unit column is sequencing id
# TODO add quast, spades, fastqc, and multiqc to yaml file
