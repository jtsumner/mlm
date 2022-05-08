# Single samples assemblies 

rule megahit:
    input:
        r1_clean = "results/bwa_out/{sample}/{sample}.fastp_bwa.r1.fastq",
        r2_clean = "results/bwa_out/{sample}/{sample}.fastp_bwa.r2.fastq"
    output:
        scaffolds = "results/megahit_out/{sample}/{sample}.contigs.fa"
    params:
        outdirec = "results/megahit_out/{sample}"
    threads: 20
    resources:
        mem="100g",
        time="10:00:00"
    shell:
        """
        module load megahit/1.0.6.1
        megahit -t {threads} -m 0.9 -1 {input.r1_clean} -2 {input.r2_clean} --out-prefix {wildcards.sample} -o {params.outdirec}_tmp
        mv {params.outdirec}_tmp/* {params.outdirec}
        rmdir {params.outdirec}_tmp
        """
