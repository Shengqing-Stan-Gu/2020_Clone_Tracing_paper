library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("AnnotationDbi")
library("org.Mm.eg.db")
library("EnsDb.Mmusculus.v79")


setwd("~/Desktop/Work Documents/Liu Lab Data analysis/2018-09_ClonTracer_RNAseq_analysis/")

rm(list=ls())

## Read in the rsem output (raw count)
rsem_count = read.csv("./RSEM_result/RSEM_rawcount.csv",row.names=1)
colnames(rsem_count) = c("14L","14R","16R","17R","20R","24L","26R","29L",
                         "30L","30R","32R","33L","33R","35L","3_invitro","4_invitro")
rsem_count[1:5,1:5]
rsem_count = apply(rsem_count, 2, as.integer)

response = c("good","good","poor","good","good","poor","poor","poor","good")

coldata = data.frame(response)
rownames(coldata) = colnames(rsem_count)[6:14]
coldata

## Convert into DESeq dataset
dds = DESeqDataSetFromMatrix(countData = rsem_count,
                                           colData = coldata,
                                           design = ~ response)



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
pdf("./Exploratory_analyses/2018-09-28_Sample_distance.pdf",height=7,width=8)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
dev.off()

## PCA plot
plotPCA(rld, intgroup=c("response"))
pcaData = plotPCA(rld, intgroup=c("response"), returnData=T)
percentVar = round(100*attr(pcaData,"percentVar"))
pdf("./Exploratory_analyses/2018-09-28_PCA.pdf",height=3.5,width=5.5)
ggplot(pcaData, aes(x=PC1, y=PC2, shape=response, color=response))+
  geom_point(size=3)+
  xlab(paste0("PC1: ", percentVar[1], "% variance"))+
  ylab(paste0("PC2: ", percentVar[2], "% variance"))
dev.off()

######################################
## Differential expression analysis ##
######################################

## Run differential expression pipeline
dds = DESeqDataSetFromMatrix(countData = rsem_count[,6:14],
                             colData = coldata,
                             design = ~ response)
dds = dds[rowSums(counts(dds))>1,]
dds = DESeq(dds)
res = results(dds, contrast=c("response","poor","good"))
res
summary(res)


## Add gene symbol information
res$symbol = mapIds(org.Mm.eg.db,
                        keys=row.names(res),
                        column="SYMBOL",
                        keytype="ENSEMBL",
                        multiVals="first")

res = res[order(res$pvalue),]



## Output comparison results
write.csv(res,"Exploratory_analyses/2018-09-28_Differential_Expression_poor_good.csv")



## Output file for GSEA
symbol = toupper(res$symbol)
log2FC = res$log2FoldChange

write.table(cbind(symbol,log2FC),"./Exploratory_analyses/2018-09-28_Gene_expr.rnk",
            sep="\t",quote=F,row.names=F)






