#### Rule to assemble concatenated sample reads with Spades ####
rule spades:
    input:
        r1_clean = "results/bwa_out/{sample}/{sample}.fastp_bwa.r1.fastq",
        r2_clean = "results/bwa_out/{sample}/{sample}.fastp_bwa.r2.fastq"
    output:
        scaffolds="results/spades_out/{sample}/scaffolds.fasta"
    params:
        out_dir=directory("results/spades_out/{sample}")

    threads: 40
    resources:
        mem="100g",
        time="10:00:00"
    shell:
        """
        module load spades/3.14.1
        spades.py -1 {input.r1_clean} -2 {input.r2_clean} -o {params.out_dir} -t {threads} -m 100 --meta
        """
