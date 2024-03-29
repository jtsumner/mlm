
library(tidyverse)
library(ggplot2)

library(ggpubr)
library(cowplot)
library(Nonpareil)

```{r}
# Get list of NonPareil outputs
np_list <- list.files(path = "../../results/nonpareil_out", pattern='*\\.npo', recursive=TRUE, full.names=TRUE)
np_length <- length(np_list)
```

# Graph NonPareil Curves
```{r}
#pdf(file="../../results/notebook_out/00_NonPareilPlot.pdf")
Nonpareil.legend(Nonpareil.curve.batch(np_list, 
                                       plot.opts = list(legend = NA, 
                                                        plot.model = F, 
                                                        plot.diversity = F,
                                                        main = NA)), cex = 0.10)
Nonpareil.curve.batch(np_list, plot.opts = list(legend = NA, 
                                                        plot.model = F, 
                                                        plot.diversity = F,
                                                        main = NA))

#dev.off()
```

# Save NonPareil Summary
```{r}
## generate summary

np_summary <-as.data.frame(
                            summary(Nonpareil.curve.batch(np_list, plot = F))
                            )
#write_csv2(np_summary, "../../results/notebook_out/00_NonPareilSummary.csv")
np_summary
```

```{r}
## compute summary stats: mean, med, sd, se
mean(summary(Nonpareil.curve.batch(np_list, plot = F))[,2])
median(summary(Nonpareil.curve.batch(np_list, plot = F))[,2])
sd(summary(Nonpareil.curve.batch(np_list, plot = F))[,2])
sd(summary(Nonpareil.curve.batch(np_list, plot = F))[,2])/sqrt(np_length)

```
```{r}
np_summary %>%
  ggplot( aes(x=C)) +
    geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8) + theme_q2r()

np_summary %>%
  ggplot( aes(x=diversity, y=C)) +
    geom_point() + theme_q2r()
```

