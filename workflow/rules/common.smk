from snakemake.utils import validate
import pandas as pd

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity

#configfile: "config/config.yaml"

#validate(config, schema="workfolschemas/config.schema.yaml")
DATASETS = ["Batch_01", "Batch_02", "Batch_03"]
READS = ["R1", "R2"]
samples = (
    pd.read_csv(config["samples"], sep="\t")
    .set_index("sample_name", drop=False)
    .sort_index()
)
samples.index.names = ["sample_name"]
#samples = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
#samples.index.names = ["sample"]
#validate(samples, schema="../schemas/samples.schema.yaml")