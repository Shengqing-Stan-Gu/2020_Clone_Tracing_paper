library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library(ensembldb)
library(EnsDb.Mmusculus.v79)

getwd()

## Read in the rsem output (raw count)
rsem_count = read.csv("./RSEM_result/RSEM_rawcount.csv",row.names=1)
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



##########################
## Exploratory analysis ##
##########################

## Filter the dataset to exclude too low expression
nrow(dds)
dds = dds[rowSums(counts(dds))>1,]
nrow(dds)

## Normalize the data with rlog
rld = rlog(dds, blind=F)

## Sample distances
sampleDists = dist(t(assay(rld)))
sampleDistMatrix = as.matrix(sampleDists)
colors = colorRampPalette(rev(brewer.pal(9,"Blues")))(255)
pdf("./Exploratory_analyses/2018-12-12_Sample_distance.pdf",height=7,width=8)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
dev.off()

## PCA plot
plotPCA(rld, intgroup=c("cell_type"))
pcaData = plotPCA(rld, intgroup=c("cell_type"), returnData=T)
percentVar = round(100*attr(pcaData,"percentVar"))

pdf("./Exploratory_analyses/2018-12-12_PCA.pdf",height=3.5,width=5.5)
ggplot(pcaData, aes(x=PC1, y=PC2, color=cell_type))+
  geom_point(size=3)+
  xlab(paste0("PC1: ", percentVar[1], "% variance"))+
  ylab(paste0("PC2: ", percentVar[2], "% variance"))
dev.off()



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

write.table(res_204, "Exploratory_analyses/2019-02-08_204.txt", quote=F, sep="\t")
write.table(res_264, "Exploratory_analyses/2019-02-08_264.txt", quote=F, sep="\t")

## Add gene symbol information
tx2gene = select(EnsDb.Mmusculus.v79,
                 keys = keys(EnsDb.Mmusculus.v79,keytype="GENEID"),
                 keytype = "GENEID",
                 columns = "GENENAME")

res_204$symbol = tx2gene[match(rownames(res_204), tx2gene[,1]),2]
res_264$symbol = tx2gene[match(rownames(res_264), tx2gene[,1]),2]
res_204 = res_204[order(res_204$pvalue),]
res_264 = res_264[order(res_264$pvalue),]

write.csv(res_204,"Exploratory_analyses/Line204_vs_parental.csv")
write.csv(res_264,"Exploratory_analyses/Line264_vs_parental.csv")


## Output gene list for GSEA input
res_204_exp = as.data.frame(res_204)
res_204_exp = res_204_exp[-which(is.na(res_204_exp$symbol)),]
res_204_exp = res_204_exp[-which(res_204_exp$symbol==""),]
res_264_exp = as.data.frame(res_264)
res_264_exp = res_264_exp[-which(is.na(res_264_exp$symbol)),]
res_264_exp = res_264_exp[-which(res_264_exp$symbol==""),]

gsea_204 = res_204_exp[,c("symbol","log2FoldChange")]
gsea_204$symbol = toupper(gsea_204$symbol)
gsea_264 = res_264_exp[,c("symbol","log2FoldChange")]
gsea_264$symbol = toupper(gsea_264$symbol)
head(gsea_204)

write.table(gsea_204,"Exploratory_analyses/2018-12-13_204_vs_P.rnk",
            quote=F, row.names=F, sep="\t")
write.table(gsea_264,"Exploratory_analyses/2018-12-13_264_vs_P.rnk",
            quote=F, row.names=F, sep="\t")


## Output gene list for David enrichment analysis
head(res_204_exp,4)
res_204_significant = res_204_exp[which(res_204_exp$pvalue<0.0297),]
head(res_264_exp,4)
res_264_significant = res_264_exp[which(res_264_exp$pvalue<0.0211),]

write.csv(res_204_significant[order(res_204_significant$log2FoldChange),],
          "Exploratory_analyses/2018-12-14_204_vs_P_significant.csv")
write.csv(res_264_significant[order(res_264_significant$log2FoldChange),],
          "Exploratory_analyses/2018-12-14_264_vs_P_significant.csv")


## Add gene symbol to TPM file
TPM_count = read.csv("./RSEM_result/RSEM_scaledTPM.csv",row.names=1)
TPM_count[1:4,1:4]
tx2gene = select(EnsDb.Mmusculus.v79,
                 keys = keys(EnsDb.Mmusculus.v79,keytype="GENEID"),
                 keytype = "GENEID",
                 columns = "GENENAME")
TPM_count$symbol = tx2gene[match(rownames(TPM_count), tx2gene[,1]),2]
TPM_count[1:4,]

write.csv(TPM_count,"./RSEM_result/2018-12-20_TPM_with_Gene_Symbol.csv")
