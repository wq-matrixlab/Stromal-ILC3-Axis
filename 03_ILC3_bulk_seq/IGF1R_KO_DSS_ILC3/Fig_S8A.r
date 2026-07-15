library(dplyr)
library(ggplot2)
library(openxlsx) 
library(stringr)

data <- read.xlsx("GSEA_analysis.xlsx")  

data$idx <- ifelse(data$NES > 0, "Up", "Down")
data$Description <- str_wrap(data$Description, width = 45)

data <- arrange(data, NES)
data$Description <- factor(data$Description, levels = data$Description)


text_gap <- max(abs(data$NES)) * 0.02 


expand_factor <- 1
y_limit <- max(abs(data$NES)) * expand_factor

p <- ggplot(data = data, aes(x = Description, y = NES, fill = idx)) + 
  geom_col(alpha = 0.8, width = 0.7) + 
  theme_bw() + 

  scale_y_continuous(
    limits = c(-y_limit, y_limit),
    breaks = pretty(c(-max(abs(data$NES)), max(abs(data$NES))), n = 5) 
  ) + 
  ylab("Normalized Enrichment Score (NES)") + 
  xlab("") + 
  
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(),
    axis.line.x = element_line(colour = "black"),
    legend.position = "top",
    legend.title = element_blank()
  ) + 

  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) + 
  coord_flip() + 
  
  scale_fill_manual(
    values = c("Up" = "#74c2d7", "Down" = "#ec8574"),
    labels = c("Up-regulated", "Down-regulated")
  ) +
  

  geom_text(
    data = filter(data, idx == "Up"),
    aes(y = -text_gap, label = Description), 
    hjust = 1, 
    size = 3.5, 
    color = "black",
    lineheight = 0.8
  ) +
    

  geom_text(
    data = filter(data, idx == "Down"),
    aes(y = text_gap, label = Description), 
    hjust = 0, 
    size = 3.5, 
    color = "black",
    lineheight = 0.8
  )


print(p)