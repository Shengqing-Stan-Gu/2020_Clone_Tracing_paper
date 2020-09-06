library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library(ensembldb)
library(EnsDb.Mmusculus.v79)

getwd()

## Read in the rsem output (raw count)
rsem_count = read.csv("ATAC/all_RP.csv",row.names=1)
colnames(rsem_count) = c("P_1","P_2","204_1","204_2","264_1","264_2")

rsem_count[1:5,1:5]
rsem_count = floor(rsem_count)
cell_type = c("P","P","204","204","264","264")

coldata = data.frame(cell_type)
rownames(coldata) = colnames(rsem_count)
coldata

## Convert into DESeq dataset
dds = DESeqDataSetFromMatrix(countData = rsem_count,
                             colData = coldata,
                             design = ~ cell_type)



######################################
## Differential expression analysis ##
######################################

## Run differential expression pipeline
dds = DESeqDataSetFromMatrix(countData = rsem_count,
                             colData = coldata,
                             design = ~ cell_type)
dds = dds[rowSums(counts(dds))>1,]
dds = DESeq(dds)
res_204 = results(dds, contrast=c("cell_type","204","P"))
summary(res_204)
res_264 = results(dds, contrast=c("cell_type","264","P"))
summary(res_264)

write.csv(res_204, "ATAC/2019-06-26_204_RP.csv")
write.csv(res_264, "ATAC/2019-06-26_264_RP.csv")


##################################
## Compare RNA-seq and ATAC-seq ##
##################################

ATAC_204 = read.csv("ATAC/2019-06-26_204_RP.csv", stringsAsFactors=F, row.names=1)
ATAC_264 = read.csv("ATAC/2019-06-26_264_RP.csv", stringsAsFactors=F, row.names=1)
RNA_204 = read.csv("Exploratory_analyses/Line204_vs_parental.csv", stringsAsFactors=F, row.names=1)
RNA_264 = read.csv("Exploratory_analyses/Line264_vs_parental.csv", stringsAsFactors=F, row.names=1)
RNA_204 = RNA_204[!duplicated(RNA_204$symbol),]
RNA_264 = RNA_264[!duplicated(RNA_264$symbol),]

common_gene = intersect(rownames(ATAC_204), RNA_204$symbol)

ATAC_204 = ATAC_204[common_gene,]
ATAC_264 = ATAC_264[common_gene,]
RNA_204 = RNA_204[which(RNA_204$symbol %in% common_gene),]
RNA_264 = RNA_264[which(RNA_264$symbol %in% common_gene),]
rownames(RNA_204) = RNA_204$symbol
rownames(RNA_264) = RNA_264$symbol
RNA_204 = RNA_204[common_gene,]
RNA_264 = RNA_264[common_gene,]

cor.test(RNA_204$log2FoldChange,ATAC_204$log2FoldChange)
cor.test(RNA_264$log2FoldChange,ATAC_264$log2FoldChange)

plot(RNA_204$log2FoldChange,ATAC_204$log2FoldChange,pch=19, cex=0.5, xlim=c(-4,4), ylim=c(-2.5,2.5))
plot(RNA_264$log2FoldChange,ATAC_264$log2FoldChange,pch=19, cex=0.5, xlim=c(-4,4), ylim=c(-2.5,2.5))

comparison_204 = data.frame(cbind(RNA_204$log2FoldChange,ATAC_204$log2FoldChange))
comparison_264 = data.frame(cbind(RNA_264$log2FoldChange,ATAC_264$log2FoldChange))
colnames(comparison_204) = c("Relative_Expression", "Relative_Regulatory_Potential")
colnames(comparison_264) = c("Relative_Expression", "Relative_Regulatory_Potential")

pdf("ATAC/2019-06-26_Correlation_204.pdf", width=5, height=4)
ggplot(comparison_204,aes(x=Relative_Expression,y=Relative_Regulatory_Potential))+
  geom_point(size=0.8)+
  xlim(-3.5,3.5)+
  ylim(-2.5,2.5)
dev.off()

pdf("ATAC/2019-06-26_Correlation_264.pdf", width=5, height=4)
ggplot(comparison_264,aes(x=Relative_Expression,y=Relative_Regulatory_Potential))+
  geom_point(size=0.8)+
  xlim(-3.5,3.5)+
  ylim(-2.5,2.5)
dev.off()


