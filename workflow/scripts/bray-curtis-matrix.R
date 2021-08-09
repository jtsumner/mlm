library("vegan")
library("tidyverse")
library("pheatmap")
data("BCI")
head("BCI")


species <- read_tsv("metaphlan/merged_abundance_table.species.allDatasets.txt")

# Transform to sample x spp. abundance matrix for vegan
species_mat <- species %>%
  gather(key = "sample_name", value = "abundance", 2:last_col()) %>%
  spread(key = "sample", value = "abundance") %>%
  column_to_rownames("sample_name") %>%
  as.matrix()

BC <- vegdist(species_mat, method = "bray") %>% as.matrix

pheatmap(BC, cluster_rows = FALSE, cluster_cols = FALSE, cellwidth = 10, cellheight = 10)

