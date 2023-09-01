
rule testr:
    input:
        "results/fastqc_out/multiqc_report.html"
    output:
        "results/mlm_out/test.txt"
    conda:
        "../envs/r.yml"
    script:
        "../scripts/script.R"