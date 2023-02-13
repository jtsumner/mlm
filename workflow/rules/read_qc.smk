rule fastqc_raw:
    input: 
        r1 = get_r1,
        r2 = get_r2
    output:
        directory("results/fastqc_out/raw_qc/{sample}")
    params:
        out_dir = "results/fastqc_out/raw_qc/{sample}"
    threads: 12
    shell:
        """
        module load fastqc/0.11.5
        fastqc -t {threads} {input.r1} {input.r2} --outdir {params.out_dir}
        """


rule fastqc_fastp:
    input: 
        r1_filtered = "results/fastp_out/{sample}/{sample}.fastp.r1.fastq.gz",
        r2_filtered = "results/fastp_out/{sample}/{sample}.fastp.r2.fastq.gz"
    output:
        "results/fastqc_out/fastp_qc/{sample}/{sample}.fastp.r1_fastqc.html",
        "results/fastqc_out/fastp_qc/{sample}/{sample}.fastp.r2_fastqc.html"
    params:
        out_dir = "results/fastqc_out/fastp_qc/{sample}"
    threads: 12
    shell:
        """
        module load fastqc/0.11.5
        fastqc -t {threads} {input.r1_filtered} {input.r2_filtered} --outdir {params.out_dir}
        """ 


rule fastqc_deconvolute:
    input: 
        r1_clean = "results/bwa_out/{sample}/{sample}.fastp_bwa.r1.fastq",
        r2_clean = "results/bwa_out/{sample}/{sample}.fastp_bwa.r2.fastq"
    output:
        "results/fastqc_out/bwa_qc/{sample}/{sample}.fastp_bwa.r1_fastqc.html",
        "results/fastqc_out/bwa_qc/{sample}/{sample}.fastp_bwa.r2_fastqc.html"
    params:
        out_dir = "results/fastqc_out/bwa_qc/{sample}"
    threads: 12
    shell:
        """
        module load fastqc/0.11.5
        fastqc -t {threads} {input.r1_clean} {input.r2_clean} --outdir {params.out_dir}
        """
