
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
        mem="80g",
        time="10:00:00"
    shell:
        """
        module load spades/3.14.1
        spades.py -1 {input.r1_clean} -2 {input.r2_clean} -o {params.out_dir} -t {threads} -m 100 --meta
        """
