### Co-assembly with megahit ###

rule concat_reads:
    input:
        cleanFastQ1 = expand("results/{dataset}/bwa/{sample}.clean.R1.fastq", zip, sample=samples["sample"], dataset=samples["dataset"]),
        cleanFastQ2 = expand("results/{dataset}/bwa/{sample}.clean.R2.fastq", zip, sample=samples["sample"], dataset=samples["dataset"])
    output:
        concatR1 = "results/allDatasets/coassembly/concat_reads/concat_reads.clean.R1.fastq",
        concatR2 = "results/allDatasets/coassembly/concat_reads/concat_reads.clean.R2.fastq"
    shell:
        """
        cat {input.cleanFastQ1} > {output.concatR1}
        cat {input.cleanFastQ2} > {output.concatR2}
        """

rule megahit_coassembly:
    input:
        concatR1 = "results/allDatasets/coassembly/concat_reads/concat_reads.clean.R1.fastq",
        concatR2 = "results/allDatasets/coassembly/concat_reads/concat_reads.clean.R2.fastq"
    output:
        scaffolds = "results/allDatasets/coassembly/megahit_result/final.contigs.fa"
    params:
        outdir = "results/allDatasets/coassembly/megahit_result/tmp"
    threads: 100
    shell:
        """
        module load megahit/1.0.6.1
        megahit -t {threads} -m 520e9 -1 {input.concatR1} -2 {input.concatR2} -o {params.outdir}
        mv {params.outdir} results/allDatasets/coassembly/
        rmdir results/allDatasets/coassembly/megahit_result
        mv results/allDatasets/coassembly/tmp results/allDatasets/coassembly/megahit_result
        """

rule quast_co:
    input:
        "results/allDatasets/coassembly/megahit_result/final.contigs.fa"
    output:
        direc=directory("results/allDatasets/coassembly/quast"),
        report="results/allDatasets/coassembly/quast/report.html",
        table=report("results/allDatasets/coassembly/quast/report.tsv", caption="report/quast_co.rst", category="ASSEMBLY", subcategory="COASSEMBLY"),
        pdf=report("results/allDatasets/coassembly/quast/report.pdf", caption="report/quast_co.rst", category="ASSEMBLY", subcategory="COASSEMBLY")

    threads: 1
    conda:
        "envs/genome_qc.yml"
    shell:
        "quast.py -o {output.direc} --threads {threads} {input}"

### For co-assembled sample assemblies 

rule index_sample_contigs:
    input:
        "results/allDatasets/single_sample_assemblies/allSamples.megahit_g1000.fa"
    output:
        "results/allDatasets/single_sample_assemblies/megahit_g1000.bwt"
    shell:
        """
        module load bwa/0.7.17
        module load samtools/1.10.1

        bwa index -p results/allDatasets/single_sample_assemblies/megahit_g1000 {input}
        """


rule map_reads_to_assembly:
    input:
        cleanFastQ1 = "results/{dataset}/bwa/{sample}.clean.R1.fastq",
        cleanFastQ2 = "results/{dataset}/bwa/{sample}.clean.R2.fastq",
        index = "results/allDatasets/single_sample_assemblies/megahit_g1000.bwt"
    output:
        sortedBam = "results/allDatasets/single_sample_assemblies/mapped_reads/{dataset}/{sample}.mapped.sorted.bam"
    params:
        genome = "results/allDatasets/single_sample_assemblies/megahit_g1000",
        sam = "results/allDatasets/single_sample_assemblies/mapped_reads/{dataset}/{sample}.mapped.sam",
        bam = "results/allDatasets/single_sample_assemblies/mapped_reads/{dataset}/{sample}.mapped.bam",
    threads: 20
    shell:
        """
        module purge all
        module load bwa/0.7.17
        module load samtools/1.10.1
        module load bedtools/2.29.2
        bwa mem -t {threads} {params.genome} {input.cleanFastQ1} {input.cleanFastQ2} > {params.sam}
        samtools view -Subh -o {params.bam} {params.sam}
        samtools sort -o {output.sortedBam} {params.bam}
        """

rule index_bam:
    input:
        sortedBam = "results/allDatasets/single_sample_assemblies/mapped_reads/{dataset}/{sample}.mapped.sorted.bam"
    output:
        bamIndex = "results/allDatasets/single_sample_assemblies/mapped_reads/{dataset}/{sample}.mapped.sorted.bam.bai"
    shell:
        """
        module load samtools/1.10.1
        samtools index {input.sortedBam}
        """


rule metabat2_depth:
    input:
        sortedBam = expand("results/allDatasets/single_sample_assemblies/mapped_reads/{dataset}/{sample}.mapped.sorted.bam",
            zip, sample=samples["sample"], dataset=samples["dataset"]),
        contigs = "results/allDatasets/single_sample_assemblies/allSamples.megahit_g1000.fa",
        bamIndex = expand("results/allDatasets/single_sample_assemblies/mapped_reads/{dataset}/{sample}.mapped.sorted.bam.bai",
            zip, sample=samples["sample"], dataset=samples["dataset"])
    output:
        depth_fi = "results/allDatasets/single_sample_assemblies/metabat2/allSamples.megahit_g1000.fa.depth.txt"
    conda:
        "envs/metabat2.yml"
    threads: 20
    shell:
        """
        jgi_summarize_bam_contig_depths \
        --outputDepth {output.depth_fi} \
        --percentIdentity 97 \
        --minContigLength 1000 \
        --minContigDepth 1.0 \
        --referenceFasta {input.contigs} {input.sortedBam}
        """

rule metabat2_bin:
    input:
        contigs = "results/allDatasets/single_sample_assemblies/allSamples.megahit_g1000.fa",
        depth_fi = "results/allDatasets/single_sample_assemblies/metabat2/allSamples.megahit_g1000.fa.depth.txt"
    output:
        bin_one = "results/allDatasets/single_sample_assemblies/metabat2/bins/bin.1.fa",
        bin_dir = directory("results/allDatasets/single_sample_assemblies/metabat2/bins/")
    conda:
        "envs/metabat2.yml"
    threads: 20
    shell:
        """
        metabat2 -t {threads} \
        --inFile {input.contigs} \
        --outFile {output.bin_dir}/bin \
        --abdFile {input.depth_fi}
        """
    

rule checkm_analysis:
    input:
        bin_dir = directory("results/allDatasets/single_sample_assemblies/metabat2/bins")
    output:
        checkm_dir = directory("results/allDatasets/single_sample_assemblies/metabat2/checkm"),
        checkm_fi = "results/allDatasets/single_sample_assemblies/metabat2/checkm/checkm_output.txt"
    threads: 12
    shell:
        """
        module load checkm/1.0.7
        checkm lineage_wf \
        --threads {threads} \
        --extension 'fa' \
        --file {output.checkm_fi} \
        {input.bin_dir} {output.checkm_dir}
        """

# to do: add config data for checkm
# add checkm plots

rule checkm_plots:
    input:
        checkm_dir = directory("results/allDatasets/single_sample_assemblies/metabat2/checkm"),
        bin_dir = directory("results/allDatasets/single_sample_assemblies/metabat2/bins")
    output:
        plots_dir = directory("results/allDatasets/single_sample_assemblies/metabat2/plots"),
        qa_plot = "results/allDatasets/single_sample_assemblies/metabat2/plots/bin_qa_plot.png"
    shell:
        """
        module load checkm/1.0.7
        checkm bin_qa_plot \
        --extension 'fa' \
        {input.checkm_dir} {input.bin_dir} {output.plots_dir}
        """

# add single sample binning
# add mash 
