library(DESeq2)
library(pheatmap)
library(ggplot2)

eg_count <-read.table("count.txt", header = T,sep="\t", comment.char = "", check.names = F)
index=order(rowMeans(eg_count[,-1]),decreasing = T)
expr_ordered=eg_count[index,]
keep=!duplicated(expr_ordered$Symbol)
expr_max=expr_ordered[keep,]
write.table (expr_max, file ="test.txt", row.names =T, col.names =T)
eg_count <-read.table("test.txt", header = T,sep="\t", row.names = 1, comment.char = "", check.names = F)

eg_countData<-eg_count[apply(eg_count,1,sum)>0,]
condition <- factor(c(rep("SL",3),rep("Ctrl",3)))
colData <- data.frame(row.names=colnames(eg_countData), condition)

dds <- DESeqDataSetFromMatrix(countData = eg_countData, colData = colData, design = ~ condition)
dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE)
res <- results(dds1,contrast = c("condition","SL","Ctrl"))
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
write.table (res1, file ="DEGresult.txt", sep = "", row.names =T, col.names =T, quote =F)