#### Rule to assemble concatenated sample reads with Spades ####
rule spades:
    input:
        r1_clean = "results/bwa_out/{sample}/{sample}.fastp_bwa.r1.fastq",
        r2_clean = "results/bwa_out/{sample}/{sample}.fastp_bwa.r2.fastq",
    output:
        direc=directory("data/spades_assemblies/{sample}"),
        scaffolds="../data/spades_assemblies/{sample}/scaffolds.fasta"
    threads: 40
    shell:
        """
        module load spades; spades.py -1 {input.r1} -2 {input.r2} -o {output.direc} -t {threads} -m 100 --meta
        """
