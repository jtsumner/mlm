
import glob
import pandas as pd
from snakemake.utils import validate


rule index_sample_contigs:
    input:
        "../results/allDatasets/single_sample_assemblies/allSamples.megahit_g1000.fa"
    output:
        "../results/allDatasets/single_sample_assemblies/megahit_g1000.bwt"
    shell:
        """
        module load bwa/0.7.17
        module load samtools/1.10.1

        bwa index -p megahit_g1000 {input}
        """

