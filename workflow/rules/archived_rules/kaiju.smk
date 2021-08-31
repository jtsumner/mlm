# For kaiju binning from bwa deconvulution pipeline

# Add the following lines to rule_all to run
        # "../resources/kaiju_head/kaiju-v1.8.0-linux-x86_64.tar.gz", # Kaiju Tar
        # "../resources/kaiju_head/kaijuDB/kaiju_db_refseq_2021-02-26.tgz", # Kaiju DB download
        # expand("../results/{dataset}/abundance/kaiju_refseq/{sample}.kaiju_refseq.txt",
        #     zip, sample=samples["sample"], dataset=samples["dataset"]), # Kaiju refseq outputs
        # "../results/allDatasets/kaiju/kaiju_abundance_table.species.allDatasets.txt", # Kaiju refseq summary


### Setup and Execute Kaiju ###
 
rule kaiju_setup:
    output:
        tar = "../resources/kaiju_head/kaiju-v1.8.0-linux-x86_64.tar.gz",
        kaijuDir = directory("../resources/kaiju_head/kaijuDir"),
    params:
        kaiju_head = "../resources/kaiju_head",
        kaiju_archive = "https://github.com/bioinformatics-centre/kaiju/releases/download/v1.8.0/kaiju-v1.8.0-linux-x86_64.tar.gz",
        kaiju_old_dir = "kaiju-v1.8.0-linux-x86_64-static"
    threads: 10
    shell:
        """
        wget {params.kaiju_archive} -P {params.kaiju_head}
        tar -xvzf {params.kaiju_head}/kaiju-v1.8.0-linux-x86_64.tar.gz -C {params.kaiju_head}
        mv {params.kaiju_head}/{params.kaiju_old_dir} {output.kaijuDir}
        """

rule kaiju_db:
    input:
        "../resources/kaiju_head/kaiju-v1.8.0-linux-x86_64.tar.gz"
    output:
        tar = "../resources/kaiju_head/kaijuDB/kaiju_db_refseq_2021-02-26.tgz",
        kaijuDB = directory("../resources/kaiju_head/kaijuDB")
    params:
        kaiju_head = "../resources/kaiju_head",
        database = "https://kaiju.binf.ku.dk/database/kaiju_db_refseq_2021-02-26.tgz"
    threads: 10
    shell:
        """
        wget {params.database} -P {output.kaijuDB}
        tar -xvzf {output.tar} -C {output.kaijuDB}
        """


rule kaiju_refseq:
    input:
        kaiju_sb = rules.kaiju_db.output.tar,
        cleanFastQ1 = "../results/{dataset}/bwa/{sample}.clean.R1.fastq",
        cleanFastQ2 = "../results/{dataset}/bwa/{sample}.clean.R2.fastq"
    output:
        profile = "../results/{dataset}/abundance/kaiju_refseq/{sample}.kaiju_refseq.txt",
    params:
        db_path = "../resources/kaiju_head/kaijuDB",
        mode = "mem"
    threads: 25
    shell:
        """
        ../resources/kaiju_head/kaijuDir/kaiju -z {threads} \
        -t {params.db_path}/nodes.dmp \
        -f {params.db_path}/kaiju_db_refseq.fmi \
        -i {input.cleanFastQ1} \
        -j {input.cleanFastQ2} \
        -a {params.mode} \
        -o {output.profile}
        """


def kaiju_merge_inputs(wildcards):
    files = expand("../results/{dataset}/abundance/kaiju_refseq/{sample}.kaiju_refseq.txt",
        zip, sample=samples["sample"], dataset=samples["dataset"])
    return files


rule kaiju_merge:
    input: 
        kaiju_merge_inputs
    output:
        "../results/allDatasets/kaiju/kaiju_abundance_table.species.allDatasets.txt"
    params:
        db_path = "../resources/kaiju_head/kaijuDB",
        taxa = "species"
    shell:
        """
        ../resources/kaiju_head/kaijuDir/kaiju2table \
        -t {params.db_path}/nodes.dmp \
        -n {params.db_path}/names.dmp \
        -r {params.taxa} \
        -o {output} \
        {input}
        """
 
