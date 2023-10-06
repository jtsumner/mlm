
rule testr:
    input:
        "results/fastqc_out/multiqc_report.html"
    output:
        "results/mlm_out/ReadNumberSummary.tsv"
    conda:
        "../envs/r.yml"
    script:
        "../scripts/script.R"