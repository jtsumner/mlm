from snakemake.utils import validate
import pandas as pd

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity

##### load config and sample sheets #####

configfile: "../config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")


samples = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
samples.index.names = ["sample"]
#validate(samples, schema="../schemas/samples.schema.yaml")