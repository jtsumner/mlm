# Multi-Level Metagenomics

Is a user-friendly, automated metagenomics pipeline the key to making your life as a bioinformatician easier? Do you still want some choice in which tools you use for each analysis step rather a rigid selection of pre-determined tools that other pipelines use? Well look no further!

Welcome to Muti-Level Metagenomics (not to be confused with Multi-Level Marketing), a flexible analysis pipeline designed to handle various metagenomics data types and analyses steps. Use MLM's flexible, low-free configuration settings to choose from multiple tools at each major step in the analysis. Yes, you heard that right: you can use the MLM pipeline to make your own pipeline. Answer the question closest to your heart. Like, is this a pipeline, or a pyramid scheme?  

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
- [Installation](#installation)
- [Setup](#setup)
- [Execution](#execution)
- [Rules](#rules)
- [Software-versions](#software-versions)
- [Configuration settings](#configuration-settings)

This is a snakemake pipeline is designed to automate common components of shotgun metagenomic data analysis. 

Briefly, reads are trimmed, deconvoluted (- human), taxnomically defined, assembled, binned, and annotated. 

Further optimization is neccessary and functional components of metagenomic analysis have yet to be integrated. 

**TODO**
* write new prep_sample_sheet helper script for updated path-based execution
* re-test megahit + spades functionality
* come up with fun name/acronym
* probably more stuff 

-------------
# Notes on snakemake 

This pipeline has been tested using snakemake and mamba installed in a single conda environment. For more information regarding snakemake, generally, please see the documentation [here](https://snakemake.readthedocs.io/en/v7.3.8/index.html)

Software versions:
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

# Installation

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

# Setup

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

# Execution

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

# Rules 

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

# Software-versions
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
# Configuration settings

So far there are six major steps in the analysis. Each of these can be turned on or off at your desire.

```
FASTQC: True
TRIM_READS: True
DECONVOLUTE: True
METAPHLAN: True
ASSEMBLE: True
METABAT2: False
```

`FASTQC` employs fastqc on raw reads and, if selected to perform the analyses which generate them, trimmed and deconvoluted reads.
`TRIM_READS` employs fastp to trim and QC raw reads.
`DECONVOLUTE` employs bwa alignment to the human reference genome to remove probable human-derived reads
`METAPHLAN` employs metaphlan3 to determine a relative abundance profile at the genus and species level
`ASSEMBLE` employs megahit or metaspades to assemble metagenomes
`METABAT2` uses the metabat algorithm to bin genomes into metagenome assembled genomes


-------------
# Citations
ASCII art generated with this [tool](http://patorjk.com/software/taag/#p=display&h=3&v=1&f=Flower%20Power&t=MULTI-LEVEL%0A%20%20%20META-%0A%20GENOMICS%0A)

```text
,---.    ,---. ___    _  .---.,---------..-./`)              .---.      .-''-. ,---.  ,---.  .-''-.   .---.      
|    \  /    .'   |  | | | ,_|\          \ .-.')             | ,_|    .'_ _   \|   /  |   |.'_ _   \  | ,_|      
|  ,  \/  ,  |   .'  | ,-./  ) `--.  ,---/ `-' \           ,-./  )   / ( ` )   |  |   |  ./ ( ` )   ,-./  )      
|  |\_   /|  .'  '_  | \  '_ '`)  |   \   `-'`"`_ _    _ _ \  '_ '`). (_ o _)  |  | _ |  . (_ o _)  \  '_ '`)    
|  _( )_/ |  '   ( \.-.|> (_)  )  :_ _:   .---.( ' )--( ' ) > (_)  )|  (_,_)___|  _( )_  |  (_,_)___|> (_)  )    
| (_ o _) |  ' (`. _` /(  .  .-'  (_I_)   |   (_{;}_)(_{;}_(  .  .-''  \   .---\ (_ o._) '  \   .---(  .  .-'    
|  (_,_)  |  | (_ (_) _)`-'`-'|__(_(=)_)  |   |(_,_)--(_,_) `-'`-'|__\  `-'    /\ (_,_) / \  `-'    /`-'`-'|___  
|  |      |  |\ /  . \ / |        (_I_)   |   |              |        \       /  \     /   \       /  |        \ 
'--'      '--' ``-'`-''  `--------'---'   '---'              `--------``'-..-'    `---`     `'-..-'   `--------` 
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
