# To run fastqc on trimmed, concat reads

rule fastqc_trimmed:
    input:
        r1="../data/concat_trimmed/{sample}.r1.fastq.gz",
        r2="../data/concat_trimmed/{sample}.r2.fastq.gz"

    output:
        "../data/concat_trimmed/{sample}.r1_fastqc.html",
        "../data/concat_trimmed/{sample}.r2_fastqc.html"
    threads: 12
    shell:
         "module load fastqc; fastqc -t {threads} {input} --outdir '../data/concat_trimmed/'"

rule multiqc_trimmed:
    input:
        expand("../data/concat_trimmed/{sample}.r1_fastqc.html", sample=samples["sample"]),
        expand("../data/concat_trimmed/{sample}.r2_fastqc.html", sample=samples["sample"])
    output:
        "multiqc_report.html"
    shell:
        "module load multiqc; multiqc ../data/concat_trimmed/"
