rule fastqc_raw:
    input: 
        r1 = get_r1,
        r2 = get_r2
    output:
        "results/fastqc_out/raw/{sample}.raw.r1_fastqc.html",
        "results/fastqc_out/raw/{sample}.raw.r2_fastqc.html"
    params:
        out_dir = "results/fastqc_out/raw/"
    threads: 12
    shell:
        "module load fastqc/0.11.5 ; fastqc -t {threads} {input} --outdir {params.out_dir}"


rule fastqc_fastp:
    input: 
        "results/fastp_out/{sample}/{sample}.fastp.r1.fastq.gz",
        "results/fastp_out/{sample}/{sample}.fastp.r2.fastq.gz"
    output:
        "results/fastqc_out/fastp_qc/{sample}.fastp.r1_fastqc.html",
        "results/fastqc_out/fastp_qc/{sample}.fastp.r2_fastqc.html"
    params:
        out_dir = "results/fastqc_out/fastp_qc/"
    threads: 12
    shell:
        "module load fastqc/0.11.5 ; fastqc -t {threads} {input} --outdir {params.out_dir}"


rule fastqc_deconvolute:
    input: 
        r1_clean = "results/bwa_out/{sample}/{sample}.fastp_bwa.r1.fastq",
        r2_clean = "results/bwa_out/{sample}/{sample}.fastp_bwa.r2.fastq"
    output:
        "results/fastqc_out/bwa_out/{sample}.fastp_bwa.r1_fastqc.html",
        "results/fastqc_out/bwa_out/{sample}.fastp_bwa.r2_fastqc.html"
    params:
        out_dir = "results/fastqc_out/bwa_out"
    threads: 12
    shell:
        "module load fastqc/0.11.5 ; fastqc -t {threads} {input} --outdir {params.out_dir}"
