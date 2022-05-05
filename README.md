# microbiome-snakemake
This is a microbiome snakemake workflow to identify and analyze viral contigs from the built miroenvironment

Note: This repository was built using https://github.com/snakemake-workflows/cookiecutter-snakemake-workflow as a temlate

-------------

## Notes on snakemake 

This pipeline has been tested using snakemake and mamba installed in a single conda environment. Software versions:
* snakemake version 7.3.8 
* mamba 0.15.3
* conda 4.10.3


Execute to create a DAG visualization of the pipeline

```
snakemake --forceall --dag | dot -Tpdf > dag.pdf
```

Execute to create a rule graph visualization 

```
snakemake --forceall --rulegraph | dot -Tpdf > dag.pdf
```


calculate_unifrac.R from 
https://github.com/biobakery/MetaPhlAn/blob/master/metaphlan/utils/calculate_unifrac.R

-------------

## Installation

1. Create a new environment to install mamba and snakemake into

```
conda env create --name mamba -c conda-forge python3=3.8 mamba
```

2. Activate the new conda environment and install snakemake using mamba

```
conda activate mamba
mamba install -c bioconda -c conda-forge snakemake=7.3.8
```

3. Confirm that snakeake has installed.

```
snakemake --version
```

## Setup

1. Move data into the data/ subdirectory in a folder named, e.g., Batch_01

2. Go into the data/ subdirectory and use the `prep_sample_sheet.sh` script to prepare a sample sheet

```
cd data/
../workflow/scripts/prep_sample_sheet.sh 
```

3. Move the output file into the config folder and confirm in config.yaml file that the samples option is correctly set to that file name

```
mv samples_minimal.tsv ../congfig/samples_minimal.tsv
```

4. Configure the cluster config file. Adjust the `account` and `partition` setting under `default-resources` to fit your cluster.

5. Open the Snakefile and adjust the rule all output to fit your desired output. 
6. Configure the general snakemake config `config/config.yaml` so that rules you want to execute are set to `True` and rules you do not want to execute are set to `False`. E.g., The following lines will assemble, perform metaphlan read based analysis, and metabat2 binning

```
ASSEMBLE: True
METAPHLAN: True
METABAT2: True
```

## Execute

1. (Optional) Check that snakemake is correctly interpretting your sample spreadsheet by executing a dryrun or one of the commands in the notes above

```
snakemake --dry-run
```

2. In the base snakemake/project directory, use the `run_snakemake.sh` file to submit the snakemake scheduler as a job to slurm 

```
sbatch run_snakeake.sh
```

alternatively, you can run an interactive and execute the contents of the `run_snakemake.sh` file if you are running into challenges with slurm

3. wait for data :)
(hopefully)

-------------

## Rules 

Generalized commands + software used by this pipeline

**FastP** 

```
fastp \
	-i {input.r1} \
	-I {input.r2} \
	--out1 {output.r1Filtered} \
	--out2 {output.r2Filtered} \
	--detect_adapter_for_pe \
	--thread {threads} \
	--length_required 50 \
	-j {output.json} \
	-h {output.html} \
	-V
```

**FastQC**

```
fastqc -t {threads} \
	{input} \
	--outdir {params.outDir}"
```

**BWA - Remove human-derived reads**

```
bwa mem -t {threads} {input.genome} {input.r1Filtered} {input.r2Filtered} > {params.sam}
samtools view -Subh -o {params.bam} {params.sam}
samtools sort -o {params.sortedBam} {params.bam}
samtools view -b -f 12 -F 256 -o {params.unmappedBam} {params.sortedBam}
bedtools bamtofastq -i {params.unmappedBam} -fq {output.cleanFastQ1} -fq2 {output.cleanFastQ2}
```

**Metaphlan**

*Setup*
```
metaphlan --install \
	--index {params.metaphlan_idx} \
	--bowtie2db {output.metaphlan_db} \
	--nproc {threads}
```

*Execute*
```
metaphlan {input.cleanFastQ1},{input.cleanFastQ2} \
        --bowtie2out {output.bowtie_out} \
        --index {params.metaphlan_idx} \
        --bowtie2db {input.metaphlan_db} \
        --nproc {threads} \
        --input_type fastq \
        --unknown_estimation \
        -o {output.profile}
```

**MegaHit**

```
megahit -t {threads} \
	-m 0.9 \
	-1 {input.cleanR1} \
	-2 {input.cleanR1} \
	-o {params.outdir_tmp}
```

**Quast**

```
quast.py \
	-o {output.direc} \
	--threads {threads} \
	--min-contig 0 \
	-L {input}
```

**Metabat2**

*Depth*
```
jgi_summarize_bam_contig_depths \
        --outputDepth {output.depth_fi} \
        --percentIdentity 97 \
        --minContigLength 1000 \
        --minContigDepth 1.0 \
        --referenceFasta {input.contigs} {input.sortedBam}
```

*Bin*
```
metabat2 -t {threads} \
        --inFile {input.contigs} \
        --outFile {output.bin_dir}/bin \
        --abdFile {input.depth_fi}
```
