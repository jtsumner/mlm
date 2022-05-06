# Single samples assemblies 

rule megahit:
    input:
        r1_clean = "results/bwa_out/{sample}/{sample}.fastp_bwa.r1.fastq",
        r2_clean = "results/bwa_out/{sample}/{sample}.fastp_bwa.r2.fastq"
    output:
        outdirec = "results/megahit_out/{sample}",
        scaffolds = "results/megahit_out/{sample}/{sample}.contigs.fa"
    threads: 20
    resources:
        mem="50g",
        time="10:00:00"
    shell:
        """
        module load megahit/1.0.6.1
        megahit -t {threads} -m 0.9 -1 {input.r1_clean} -2 {input.r2_clean} --out-prefix {wildcards.sample} -o {output.outdirec}_tmp
        mv {output.outdirec}_tmp {output.outdirec}

        """
