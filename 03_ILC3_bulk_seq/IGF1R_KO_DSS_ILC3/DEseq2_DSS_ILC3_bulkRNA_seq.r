setwd("D:/IGF1_bioinfo/bulkRNA_DSS_ILC3")

library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(clusterProfiler)
library(openxlsx)

data = read.csv("count.csv", header = TRUE, sep = ",")

index=order(rowMeans(data[,-1]),decreasing = T)
expr_ordered=data[index,]
keep=!duplicated(expr_ordered$symbol)
expr_max=expr_ordered[keep,]

write.xlsx (expr_max, file ="expr_max.xlsx", rown.ames =T, coln.ames =T)

eg_count <-read.table("expr_max.txt", header = T,sep="\t", row.names = 1, comment.char = "", check.names = F)

eg_countData<-eg_count[apply(eg_count,1,sum)>10,]
write.xlsx (eg_countData, file ="clean_count.xlsx", row.names =T, col.names =T)

condition <- factor(c(rep("KO",4),rep("WT",4)))
colData <- data.frame(row.names=colnames(eg_countData), condition)

dds <- DESeqDataSetFromMatrix(countData = eg_countData, colData = colData, design = ~ condition)
dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE)
res <- results(dds1,contrast = c("condition","KO","WT"))
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
write.xlsx (res1, file ="DEGresult.xlsx",row.names =T, col.names =T, quote =F)