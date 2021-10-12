# microbiome-snakemake
This is a microbiome snakemake workflow to identify and analyze viral contigs from the built miroenvironment

Note: This repository was built using https://github.com/snakemake-workflows/cookiecutter-snakemake-workflow as a temlate

-------------

## Notes on snakemake 

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
