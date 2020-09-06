library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library(ensembldb)
library(EnsDb.Mmusculus.v79)
library("fgsea")

getwd()

## Read in the rsem output (raw count)
rsem_count = read.table("GSE83848_RNAseq_data_matrix_all.txt", row.names=1, header=T, stringsAsFactors=F)
rsem_count = rsem_count[,1:6]

rsem_count[1:5,1:5]
rsem_count = floor(rsem_count)
cell_type = c("NT","NT","NT","IFNg","IFNg","IFNg")

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
IFNg_sig = results(dds, contrast=c("cell_type","IFNg","NT"))
summary(IFNg_sig)

write.csv(IFNg_sig, "2019-06-25_Minn_IFNg_deseq.csv")



##############################
## Comparison with Line 264 ##
##############################

Diff_Minn = read.csv("2019-06-25_Minn_IFNg_deseq.csv",row.names=1)
Diff_264 = read.csv("Line264_vs_parental.csv",row.names=1)

common_genes = intersect(rownames(Diff_Minn), rownames(Diff_264))

head(Diff_Minn,10)
head(Diff_264,10)


Minn_sig_genes = rownames(Diff_Minn)[which(Diff_Minn$padj<0.1)]

pdf("2019-06-26_Minn_significant_genes_dotplot.pdf",width=8,height=7)
plot(Diff_Minn[Minn_sig_genes,]$log2FoldChange,Diff_264[Minn_sig_genes,]$log2FoldChange, pch=19, cex=0.5,
     xlab="Minn_log2FC", ylab="B64_log2FC")
dev.off()
cor.test(Diff_Minn[Minn_sig_genes,]$log2FoldChange,Diff_264[Minn_sig_genes,]$log2FoldChange)

pdf("2019-06-26_Minn_all_genes_dotplot.pdf",width=8,height=7)
plot(Diff_Minn[common_genes,]$log2FoldChange,Diff_264[common_genes,]$log2FoldChange, pch=19, cex=0.5,
     xlab="Minn_log2FC", ylab="B64_log2FC")
dev.off()
cor.test(Diff_Minn[common_genes,]$log2FoldChange,Diff_264[common_genes,]$log2FoldChange)


## GSEA analysis on the top and bottom z score genes regarding their expression in 204 and 264
data(examplePathways)
data(exampleRanks)
head(names(examplePathways),4)
head(exampleRanks,4)

Up_set = row.names(Diff_Minn[which(Diff_Minn$padj<0.1 & Diff_Minn$log2FoldChange>0),])
Down_set = row.names(Diff_Minn[which(Diff_Minn$padj<0.1 & Diff_Minn$log2FoldChange<0),])
Minn_sets = list(Up_set, Down_set)
names(Minn_sets) = c("Long-term IFNg-induced genes","Long-term IFNg-suppressed genes")

Rank_264 = Diff_264$log2FoldChange
names(Rank_264) = rownames(Diff_264)


GSEA_264 = fgsea(pathways = Minn_sets, 
                 stats = Rank_264,
                 minSize=5,
                 maxSize=500,
                 nperm=10000)

write.csv(GSEA_264[order(pval),1:7],"2019-06-26_GSEA_264_Minn.csv")


pdf("2019-06-26_GSEA_Minn_up.pdf",width=5,height=3)
plotEnrichment(Minn_sets[[GSEA_264[order(pval), ]$pathway[1]]],Rank_264)+
  labs(title=GSEA_264[order(pval), ]$pathway[1])
dev.off()

pdf("2019-06-26_GSEA_Minn_down.pdf",width=5,height=3)
plotEnrichment(Minn_sets[[GSEA_264[order(pval), ]$pathway[2]]],Rank_264)+
  labs(title=GSEA_264[order(pval), ]$pathway[2])
dev.off()

write.csv(Diff_264[Up_set,],"2019-08-12_Minn_up_genes_in_264.csv")
write.csv(Diff_264[Down_set,],"2019-08-12_Minn_down_genes_in_264.csv")
