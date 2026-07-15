library(ggplot2)
library(ggrepel)
deg.data <- read.table("DEG.txt",
sep = "\t",
header = TRUE,
check.names = FALSE)
cut_off_Pvalue = 0.05
cut_off_log2FC = 0.26
deg.data$change = ifelse(
deg.data$pvalue < cut_off_Pvalue & abs(deg.data$log2FoldChange) >= cut_off_log2FC,
ifelse(deg.data$log2FoldChange > cut_off_log2FC, "Up", "Down"),
"Stable"
)
p <- ggplot(
deg.data,
aes(x = log2FoldChange,
y = -log10(pvalue),
colour = change)
) +
geom_point(alpha = 0.4, size = 3.5) +
scale_color_manual(values = c(
 "Up" = "#ca0020",
  "Stable" = "gray80",
  "Down" = "#0571b0"
)) +
geom_vline(xintercept = c(-cut_off_log2FC, cut_off_log2FC),
lty = 4, col = "black", lwd = 0.8) +
geom_hline(yintercept = -log10(cut_off_Pvalue),
lty = 4, col = "black", lwd = 0.8) +
labs(x = "log2(fold change)",
y = "-log10(p-value)") +
theme_bw() +
theme(
plot.title = element_text(hjust = 0.5),
legend.position = "right",
legend.title = element_blank()
)
p + xlim(-7, 7)+ coord_cartesian(ylim = c(0, 18))
q <- p +
  geom_text_repel(
    data = deg.data,
    aes(
      x = log2FoldChange,
      y = -log10(pvalue),
      label = Label
    ),
    size = 5,
    box.padding = unit(0.9, "lines"),
    point.padding = unit(0.8, "lines"),
    segment.color = "black",
    max.overlaps = 30000,
    show.legend = FALSE
  )

q + xlim(-7, 7)+ coord_cartesian(ylim = c(0, 18))