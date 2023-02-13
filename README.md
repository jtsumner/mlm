# Multi-Level Metagenomics

Is a user-friendly, automated metagenomics pipeline the key to making your life as a bioinformatician easier? Do you still want some choice in which tools you use for each analysis step rather a rigid selection of pre-determined tools that other pipelines use? Well look no further!

Welcome to Muti-Level Metagenomics (not to be confused with Multi-Level Marketing), a flexible snakemake workflow designed to handle various metagenomics data types and analyses steps. Use MLM's flexible, high-level configuration settings to choose from multiple tools at each major step in the analysis. Yes, you heard that right: you can use the MLM pipeline to make your own pipeline that answers the questions closest to your heart. Like, is this a pipeline, or a pyramid scheme?  

**UNDER ACTIVE DEVELOPMENT** :) Use at you own risk...

```text
,---.    ,---. .---.    ,---.    ,---. 
|    \  /    | | ,_|    |    \  /    | 
|  ,  \/  ,  ,-./  )    |  ,  \/  ,  | 
|  |\_   /|  \  '_ '`)  |  |\_   /|  | 
|  _( )_/ |  |> (_)  )  |  _( )_/ |  | 
| (_ o _) |  (  .  .-'  | (_ o _) |  | 
|  (_,_)  |  |`-'`-'|___|  (_,_)  |  | 
|  |      |  | |        |  |      |  | 
'--'      '--' `--------'--'      '--' 
```

**Table of Contents**
- [Multi-Level Metagenomics](#multi-level-metagenomics)
- [Notes on snakemake](#notes-on-snakemake)
- [Getting Started](#getting-started)
	- [Installation](#installation)
		- [STEP 1: Clone Repository](#step-1-clone-repository)
		- [STEP 2: (OPTION 1) Conda Install from .yaml file](#step-2-option-1-conda-install-from-yaml-file)
		- [STEP 2: (OPTION 2) Manual Conda Install](#step-2-option-2-manual-conda-install)
	- [Setup](#setup)
	- [Execution](#execution)
- [Pipeline Contents](#pipeline-contents)
	- [Rules](#rules)
	- [Software-versions](#software-versions)
	- [Configuration settings](#configuration-settings)
- [References](#references)

This is a snakemake pipeline is designed to automate common components of shotgun metagenomic data analysis. Briefly, reads are trimmed, deconvoluted (- human), taxnomically defined, assembled, binned, and annotated. Further optimization is neccessary and functional components of metagenomic analysis have yet to be integrated. 

**TODO**
* Update prep_sample_sheet.py so it sorts samples into alphanumeric order
* check metabat
* integrate visualizations + comparisons
* move logs to results folder
* clean up bowtie2 code to reduce used disk space + gzip files
* add bioconda biopython=1.78  to snakemamba base install
* probably more stuff 

**MaybeTODO**
* full setup/install script
* singularity
* test cases?

-------------
# Notes on snakemake 

This pipeline has been tested using snakemake and mamba installed in a single conda environment. For more information regarding snakemake, generally, please see the documentation [here](https://snakemake.readthedocs.io/en/v7.3.8/index.html)

Execute to create a DAG visualization of the pipeline

```
snakemake --forceall --dag | dot -Tpdf > dag.pdf
```

Execute to create a rule graph visualization 

```
snakemake --forceall --rulegraph | dot -Tpdf > dag.pdf
```

-------------
# Getting Started 

## Installation

Tested versions for base conda software environment:
* snakemake version 7.3.8 
* mamba 0.15.3
* anaconda 4.10.3
### STEP 1: Clone Repository
1. Clone the MLM repository and move into the freshly cloned directory
	
	```
	git clone https://github.com/jtsumner/mlm.git
	cd mlm
	```

### STEP 2: (OPTION 1) Conda Install from .yaml file
1. (Optional) if you are on an HPC (i.e., quest), load the anaconda or miniconda module
   
	```
	module load mamba
	```

2. Create a new conda environment based on the mamba yaml file located in `workflow/envs/mamba.yml`

	```
	conda env create -f workflow/envs/mamba.yml
	```

3. The base environment with snakemake and mamba should now be available using `source activate mamba`

### STEP 2: (OPTION 2) Manual Conda Install

1. Create a new conda environment to install mamba and snakemake. If you have not already installed conda, install it using the documentatio found [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

	```
	conda create --name mamba -c conda-forge -c bioconda
	```

2. Activate the new conda environment and install snakemake using mamba

	```
	mamba create -c conda-forge -c bioconda -n snakemake snakemake
	```

3. Confirm that snakeake has installed.

	```
	snakemake --version
	```

## Setup

1. Move data into the `data/` subdirectory

2. From the mlm base directory (mlm or metagenomics-snsakemake), use the `prep_sample_sheet.py` script in the `workflow/scripts` subdirectory to prepare a sample sheet. It will automatically create a file called `sample_sheet.tsv` in the `config/` subdirect E.g.,

	```
	python3 workflow/scripts/prep_sample_sheet.py --delimiter="_S"
	```

3. Configure the **cluster config** file. Adjust the `account` and `partition` setting under `default-resources` to fit your cluster. Note that the current cluster configurations is a basic setup for SLURM on Quest that's based on this [blog](https://fame.flinders.edu.au/blog/2021/08/02/snakemake-profiles-updated)

4. Open the Snakefile and adjust the rule all output to fit your desired output. 
   
5. [Configure](#configuration-settings) the general snakemake config `config/config.yaml` so that rules you want to execute are set to `True` and rules you do not want to execute are set to `False`. E.g., The following lines will assembl and perform metaphlan read-based analysis, but won't execute metabat2 binning. 

	```
	FASTQC: False
	NONPAREIL: False
	TRIM_READS: True
	DECONVOLUTE: True
	BOWTIE2: True
	METAPHLAN: False
	ASSEMBLE: False
	MEGAHIT: True
	SPADES: True
	METABAT2: False
	```
6. Download Bowtie2 index for human reference genome
```
mkdir resources/bowtie_human
cd resources/bowtie_human
wget https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip
unzip GRCh38_noalt_as.zip
```

## Execution

1. (Optional) Check that snakemake is correctly interpretting your sample spreadsheet by executing a dryrun or one of the commands in the notes above

	```
	snakemake --dry-run
	```

2. In the base snakemake/project directory, use the `run_snakemake.sh` file to submit the snakemake scheduler as a job to slurm. (Alternatively, you can run an interactive job snakeake via command line for troubleshooting)

	```
	sbatch run_snakeake.sh
	```

3. wait for data :)
(hopefully)

### Use interactive jobs on SLURM for debugging
Start interactive job on slurm

```
srun -A b1042 --partition=genomics -N 1 -n 24 --mem=64G --time=02:00:00 --pty bash -i
```

Start interactive job on slurm with salloc 
ssh onto the qnode it outputs
```
salloc -A p31648 --partition=long -N 1 -n 24 --mem=64G --time=12:00:00 bash -i
```

-------------
# Pipeline Contents

## Rules 

Generalized commands + software used by this pipeline

**FastP** 

```
fastp -i {input.r1} \
        -I {input.r2} \
        --out1 {output.r1_filtered} \
        --out2 {output.r2_filtered} \
        --detect_adapter_for_pe \
        --dedup \
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

**Bowtie2 - Remove human-derived reads**
```
        bowtie2 -p {threads} -x {params.filter_db} --very-sensitive -1 {input.r1} -2 {input.r2}| \
        samtools view -bS -@ {threads}| \
        samtools sort -@ {threads} -n -o {output.sorted_bam}

        samtools fastq -1 {output.r1_clean} -2 {output.r2_clean} -@ {threads} -f 12 -F 256 {output.sorted_bam}

        samtools flagstat -@ {threads} -O tsv {output.sorted_bam} > {output.flagstat}
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

## Software-versions
* bedtools v2.29.2
* biopython v1.78
* bowtie2 v2.4.4
* bwa v0.7.17
* checkm v1.0.7
* fastp v0.20.1
* fastqc v0.11.9
* fastqc v0.11.5
* hclust2 v1.0.0
* kneaddata v0.10.00
* megahit v1.0.6.1
* metabat2 v2.15
* metaphlan v3.0.13
* python v2.7
* python v3.7
* quast v5.0.2
* samtools v1.10.1
* spades v3.14.1

to add and/or deprecated:
* kaiju v1.8.0
* vcontact2
* blast
* diamond
* mcl
* hmmer
* cython v0.29.21
* scikit-learn v0.21.3* prodigal

-------------
## Configuration settings

So far there are six major steps in the analysis. Each of these can be turned on or off at your desire.

```
FASTQC: True
TRIM_READS: True
DECONVOLUTE: True
METAPHLAN: True
ASSEMBLE: True
METABAT2: False
```

* `FASTQC` employs fastqc on raw reads and, if selected to perform the analyses which generate them, trimmed and deconvoluted reads.
* `TRIM_READS` employs fastp to trim and QC raw reads.
* `DECONVOLUTE` employs bwa alignment to the human reference genome to remove probable human-derived reads
* `METAPHLAN` employs metaphlan3 to determine a relative abundance profile at the genus and species level
* `ASSEMBLE` employs megahit or metaspades to assemble metagenomes
* `METABAT2` uses the metabat algorithm to bin genomes into metagenome assembled genomes


-------------
# References

* Calculate_unifrac.R from metaphlan/biobakery [here](https://github.com/biobakery/MetaPhlAn/blob/master/metaphlan/utils/calculate_unifrac.R)
* ASCII art generated with this [tool](http://patorjk.com/software/taag/#p=display&h=3&v=1&f=Flower%20Power&t=MULTI-%0ALEVEL%0AMETA-%0AGENOMICS%0A)

```text
,---.    ,---. ___    _  .---.,---------..-./`)                                             
|    \  /    .'   |  | | | ,_|\          \ .-.')                                            
|  ,  \/  ,  |   .'  | ,-./  ) `--.  ,---/ `-' \                                            
|  |\_   /|  .'  '_  | \  '_ '`)  |   \   `-'`"`_ _    _ _                                  
|  _( )_/ |  '   ( \.-.|> (_)  )  :_ _:   .---.( ' )--( ' )                                 
| (_ o _) |  ' (`. _` /(  .  .-'  (_I_)   |   (_{;}_)(_{;}_)                                
|  (_,_)  |  | (_ (_) _)`-'`-'|__(_(=)_)  |   |(_,_)--(_,_)                                 
|  |      |  |\ /  . \ / |        (_I_)   |   |                                             
'--'      '--' ``-'`-''  `--------'---'   '---'                                             
  .---.      .-''-. ,---.  ,---.  .-''-.   .---.                                            
  | ,_|    .'_ _   \|   /  |   |.'_ _   \  | ,_|                                            
,-./  )   / ( ` )   |  |   |  ./ ( ` )   ,-./  )                                            
\  '_ '`). (_ o _)  |  | _ |  . (_ o _)  \  '_ '`)                                          
 > (_)  )|  (_,_)___|  _( )_  |  (_,_)___|> (_)  )                                          
(  .  .-''  \   .---\ (_ o._) '  \   .---(  .  .-'                                          
 `-'`-'|__\  `-'    /\ (_,_) / \  `-'    /`-'`-'|___                                        
  |        \       /  \     /   \       /  |        \                                       
  `--------``'-..-'    `---`     `'-..-'   `--------`                                       
,---.    ,---.   .-''-.,---------.   ____                                                   
|    \  /    | .'_ _   \          \.'  __ `.                                                
|  ,  \/  ,  |/ ( ` )   `--.  ,---/   '  \  \                                               
|  |\_   /|  . (_ o _)  |  |   \  |___|  /  | _ _    _ _                                    
|  _( )_/ |  |  (_,_)___|  :_ _:     _.-`   |( ' )--( ' )                                   
| (_ o _) |  '  \   .---.  (_I_)  .'   _    (_{;}_)(_{;}_)                                  
|  (_,_)  |  |\  `-'    / (_(=)_) |  _( )_  |(_,_)--(_,_)                                   
|  |      |  | \       /   (_I_)  \ (_ o _) /                                               
'--'      '--'  `'-..-'    '---'   '.(_,_).'                                                
  .-_'''-.      .-''-. ,---.   .--.   ,-----.   ,---.    ,---.-./`)    _______     .-'''-.  
 '_( )_   \   .'_ _   \|    \  |  | .'  .-,  '. |    \  /    \ .-.')  /   __  \   / _     \ 
|(_ o _)|  ' / ( ` )   |  ,  \ |  |/ ,-.|  \ _ \|  ,  \/  ,  / `-' \ | ,_/  \__) (`' )/`--' 
. (_,_)/___|. (_ o _)  |  |\_ \|  ;  \  '_ /  | |  |\_   /|  |`-'`",-./  )      (_ o _).    
|  |  .-----|  (_,_)___|  _( )_\  |  _`,/ \ _/  |  _( )_/ |  |.---.\  '_ '`)     (_,_). '.  
'  \  '-   .'  \   .---| (_ o _)  : (  '\_/ \   | (_ o _) |  ||   | > (_)  )  __.---.  \  : 
 \  `-'`   | \  `-'    |  (_,_)\  |\ `"/  \  ) /|  (_,_)  |  ||   |(  .  .-'_/  \    `-'  | 
  \        /  \       /|  |    |  | '. \_/``".' |  |      |  ||   | `-'`-'     / \       /  
   `'-...-'    `'-..-' '--'    '--'   '-----'   '--'      '--''---'   `._____.'   `-...-'   
                                                                                            
```
