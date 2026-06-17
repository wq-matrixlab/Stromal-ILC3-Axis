setwd("D:/")

library(ggplot2)
library(ggrepel)

deg.data <- read.table("V1.txt",header=TRUE,check.names =F, sep = "\t")

cut_off_negLogPvalue = 3
cut_off_log2FC = 0.58

deg.data$change = ifelse(deg.data$negLogPvalue > cut_off_negLogPvalue & abs(deg.data$log2FC) >= cut_off_log2FC, 
                     ifelse(deg.data$log2FC> cut_off_log2FC ,'Up','Down'),
                     'Normal')

p<-ggplot(
  deg.data, 
  aes(x = log2FC, 
      y = negLogPvalue, 
      colour=change)) +
      geom_point(alpha=1, size=2) +
      scale_color_manual(values=c("#3C3CEF", "#d2dae2","#e94234"))+
  geom_vline(xintercept=c(-0.58,0.58),lty=3,col="black",lwd=1) +
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=1) +
  labs(x="log2(fold change)",
       y="negLogPvalue")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(hjust = 0.8), 
        legend.position="right", 
        legend.title = element_blank())+
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1))+
  theme(axis.text.y = element_text(size = 12, angle = 45, hjust = 1))+
  theme(axis.title=element_text(size=12))
  
scale_x_continuous(
    limits = c(-9, 9),
    breaks = seq(-9, 9, by = 3)
  )
p
q<-p+geom_text_repel(data = deg.data, aes(x = log2FC, 
                                      y = negLogPvalue, 
                                      label = Label),
                  size = 5,box.padding = unit(1.01, "lines"),
                  point.padding = unit(0.8, "lines"), 
                  segment.color = "black",
                  max.overlaps =30000,				  
                  show.legend = FALSE,
				  color="black")
				  
				
q