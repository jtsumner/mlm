### Remove contaminant reads aligning to human reference genome ###
rule get_human_genome:
    output:
        "resources/genome/{}".format(config["genome_name"])
    params:
        human_genome = config["human_genome"]
    resources:
        time="02:00:00"
    shell:
        """
        cd resources/genome/
        wget {params.human_genome}
        """

rule index_human_genome:
    input:
        "resources/genome/{}".format(config["genome_name"])
    output:
        "resources/genome/{}.ann".format(config["genome_name"])
    params:
        genome_index = config["genome_name"]
    threads: 10
    resources:
        mem="10G",
        time="04:00:00"
    shell:
        """
        module load bwa/0.7.17
	    bwa index -a bwtsw {input} {params.genome_index}
        """

rule bwa_map:
    input:
        r1_filtered = "results/fastp_out/{sample}/{sample}.fastp.r1.fastq.gz",
        r2_filtered = "results/fastp_out/{sample}/{sample}.fastp.r2.fastq.gz",
        genome_idxed = "resources/genome/{}.ann".format(config["genome_name"]),
        genome = "resources/genome/{}".format(config["genome_name"])
    output:
        r1_clean = "results/bwa_out/{sample}/{sample}.fastp_bwa.r1.fastq",
        r2_clean = "results/bwa_out/{sample}/{sample}.fastp_bwa.r2.fastq"
    params:
        human_genome = config["genome_name"],
        #genome = "resources/genome/{params.human_g}.ann",
        sam = "results/bwa_out/{sample}/{sample}.mapped.sam",
        bam = "results/bwa_out/{sample}/{sample}.mapped.bam",
        sorted_bam = "results/bwa_out/{sample}/{sample}.mapped.sorted.bam",
        unmapped_bam = "results/bwa_out/{sample}/{sample}.unmapped.bam"
    threads: 20
    resources:
        mem="15G",
        time="08:00:00"
    shell:
        """
        module purge all
        module load bwa/0.7.17
        module load samtools/1.10.1
        module load bedtools/2.29.2
        bwa mem -t {threads} {input.genome} {input.r1_filtered} {input.r2_filtered} > {params.sam}
        samtools view -Subh -o {params.bam} {params.sam}
        samtools sort -o {params.sorted_bam} {params.bam}

        samtools view -b -f 12 -F 256 -o {params.unmapped_bam} {params.sorted_bam}
        bedtools bamtofastq -i {params.unmapped_bam} -fq {output.r1_clean} -fq2 {output.r2_clean}
        """







#samtools view -bS -@ 20 B16_LyPMA.mapped.sam > B16_LyPMA.mappedtest.bam
#samtools view -b -@ 20 -f 12 -F 256 B16_LyPMA.mappedtest.bam > B16_LyPMA.unmappedtest.bam
#samtools sort -n -@ 20 B16_LyPMA.unmappedtest.bam -o B16_LyPMA.unmappedtest.sorted.bam
#samtools fastq -@ 20 B16_LyPMA.unmappedtest.sorted.bam -1 B16_LyPMA.unmappedtest.r1.fastq -2 B16_LyPMA.unmappedtest.r2.fastq
