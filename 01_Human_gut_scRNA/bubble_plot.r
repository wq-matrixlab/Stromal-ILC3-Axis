library(ggplot2)
library(dplyr)

go_results <- read.csv("pathway.csv")

go_results <- go_results %>%
  mutate(log_pvalue = -log10(pvalue))

ggplot(go_results, aes(x = GeneRatio, y = reorder(Description, GeneRatio))) +
    geom_point(aes(size = Count, color = log_pvalue)) +
    scale_color_gradient(low = "#8A2BE2", high = "red", limits = c(1, 10)) + 
    scale_size(range = c(2, 10)) +
    theme_bw() +
    labs(
        x = "Gene Ratio",
        y = "GO Terms",
        color = "-log10(pvalue)",
        size = "Count"
    ) +
    theme(
        axis.text.y = element_text(size = 10, angle = 0, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    ) +
    coord_cartesian(xlim = c(0.014, 0.036))  
