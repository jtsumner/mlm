# do_something <- function(data_path, out_path, threads, myparam) {
#     # R code
# }

# do_something(snakemake@input[[1]], snakemake@output[[1]], snakemake@threads, snakemake@config[["myparam"]])
source("workflow/r/nature_theme.r")
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggridges)

# Set multiqc file name and summary output name
readqc = "results/fastqc_out/multiqc_data/multiqc_fastqc.txt"
summary_out = "results/mlm_out/ReadNumberSummary.tsv"


summary <- read_tsv(readqc, name_repair = "universal") %>% 
        separate(Sample, into = c("Tool", "Sample", "Read"), sep = "\\s*\\|\\s*") %>%     
        mutate(Read = case_when(
                str_detect(Filename, "R1|r1") ~ "R1", 
                str_detect(Filename, "R2|r2") ~ "R2",  
                TRUE ~ "RM")) %>%
        mutate(ToolLevel = case_when(
                Tool == "raw_qc" ~ 1,
                Tool == "fastp_qc" ~ 2,
                Tool == "bbduk_qc" ~ 3,
                Tool == "bowtie_qc" ~ 4,
                Tool == "bbmerge_qc" ~ 5,
                TRUE ~ 6))


FoldChanges <- summary %>% select(Sample, Read, Tool, ToolLevel, Total.Sequences)

FoldChanges <- FoldChanges %>% left_join(FoldChanges, by = c("Sample", "Read"))%>%
  #filter(ToolLevel.x >= ToolLevel.y) %>%
  mutate(FoldChange = Total.Sequences.x / Total.Sequences.y,
         RelativeTool = Tool.y,
         Tool = Tool.x,
         ToolLevel = ToolLevel.x) %>%
  select(Sample, Read, Tool, ToolLevel,FoldChange, RelativeTool)


summary <- FoldChanges %>% 
  pivot_wider(names_from = RelativeTool, 
              values_from = FoldChange
              ) %>% 
  left_join(summary)


write.csv(summary, summary_out)


read_list <- unique(FoldChanges$Read)
tool_list <- sort(unique(FoldChanges$RelativeTool))

pdf("text.pdf", height=4, width=4)

for (tool in tool_list) {
  for (read in read_list) {

    print(
      FoldChanges %>%
        filter(Read == read,
              RelativeTool == tool
              ) %>%
        ggplot(., aes(x=FoldChange, y=Tool, fill=Tool)) + 
          geom_density_ridges(
            aes(point_color = Tool, 
                point_fill = Tool, 
                point_shape = Tool), 
            alpha = .2, 
            point_alpha = 1, 
            jittered_points = TRUE, 
            scale=.9, 
            size=1, 
            quantile_lines = TRUE, 
            quantiles = 2
            ) +
          nature_theme("", "") +
          labs(
            x=paste0(read, 
                    "(# of Reads) / ",
                    read, 
                    "(# of Reads at ", tool, 
                    ")"
                    ), 
            y="Processing Step"
            ) +
          theme(
            aspect.ratio = 1, 
            legend.position = "none",
            axis.ticks.y = element_blank(),
            ) + 
          scale_point_color_hue(l = 40) +
          scale_discrete_manual(
            aesthetics = "point_shape", 
            values = c(21, 22, 23, 24)
            ) +
          ggtitle(paste("Fold Chage in Num. Reads", read, tool))              
    )
  }
}

dev.off()
