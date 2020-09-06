rm(list=ls())
library("EnsDb.Mmusculus.v79")
library("EnsDb.Hsapiens.v79")
library("biomaRt")
library("fgsea")
library("ggplot2")

setwd("/Users/Shengqing/Desktop/Work Documents/Liu Lab Data analysis/2018-12-11_Clone_204_264/")

## Read in the DESeq2 result on 204 vs parental
expr_diff_204 = read.csv("Exploratory_analyses/Line204_vs_parental.csv",row.names=1)
head(expr_diff_204)
expr_diff_264 = read.csv("Exploratory_analyses/Line264_vs_parental.csv",row.names=1)
head(expr_diff_264)


## Read in the reference z scores
reference_scores = read.csv("TIDE_dysfunction/2018-12-17_Reference_z_scores.csv",row.names=1)
reference_scores$symbol = mapIds(EnsDb.Hsapiens.v79,
                                 keys=row.names(reference_scores),
                                 column="GENENAME",
                                 keytype="ENTREZID",
                                 multiVals="first")

reference_scores = reference_scores[order(reference_scores[,1]),]
reference_scores = reference_scores[-which(is.na(reference_scores$symbol)),]
head(reference_scores)

reference_scores = reference_scores[reference_scores$symbol %in% toupper(expr_diff_204$symbol),]


## Retrieve the differential expression of high-z genes and low-z genes
high_z_genes = tail(reference_scores$symbol,200)
low_z_genes = head(reference_scores$symbol,200)

## Convert human gene symbol into mouse gene symbol
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

high_z_mouse = getLDS(attributes = "hgnc_symbol",
                        filters = "hgnc_symbol",
                        values = high_z_genes,
                        mart = human,
                        attributesL = "mgi_symbol",
                        martL = mouse,
                        uniqueRows = T)
low_z_mouse = getLDS(attributes = "hgnc_symbol",
                     filters = "hgnc_symbol",
                     values = low_z_genes,
                     mart = human,
                     attributesL = "mgi_symbol",
                     martL = mouse,
                     uniqueRows = T)
human_to_mouse = getLDS(attributes = "hgnc_symbol",
                        filters = "hgnc_symbol",
                        values = reference_scores$symbol,
                        mart = human,
                        attributesL = "mgi_symbol",
                        martL = mouse,
                        uniqueRows = T)

## Compare the global association between differential expression and z score
Comparison_204 = human_to_mouse
logFC = vector()
zscore = vector()
for (i in 1:dim(Comparison_204)[1]) {
  logFC = c(logFC, expr_diff_204[which(expr_diff_204$symbol==Comparison_204$MGI.symbol[i]),2][1])
  zscore = c(zscore, reference_scores[which(reference_scores$symbol==Comparison_204$HGNC.symbol[i]),1][1])
}
Comparison_204 = cbind(Comparison_204,logFC,zscore)
cor_204 = cor.test(Comparison_204$logFC,Comparison_204$zscore, method="spearman",na.rm=T)


Comparison_264 = human_to_mouse
logFC = vector()
zscore = vector()
for (i in 1:dim(Comparison_264)[1]) {
  logFC = c(logFC, expr_diff_264[which(expr_diff_264$symbol==Comparison_264$MGI.symbol[i]),2][1])
  zscore = c(zscore, reference_scores[which(reference_scores$symbol==Comparison_264$HGNC.symbol[i]),1][1])
}
Comparison_264 = cbind(Comparison_264,logFC,zscore)
cor_264 = cor.test(Comparison_264$logFC,Comparison_264$zscore, method="spearman")

write.csv(cbind(c(cor_204$estimate,cor_204$p.value),c(cor_264$estimate,cor_264$p.value)),
            "TIDE_dysfunction/2019-06-13_B04_logFC_zscore_correlation.csv")


pdf("TIDE_dysfunction/2019-06-13_B04_logFC_zscore_dotplot.pdf",width=5,height=4)
plot(Comparison_204$logFC,Comparison_204$zscore,pch=19,cex=0.2,
     xlab="B04_RNAseq_logFC",ylab="TIDE_z_score")
dev.off()

pdf("TIDE_dysfunction/2019-06-13_B64_logFC_zscore_dotplot.pdf",width=5,height=4)
plot(Comparison_264$logFC,Comparison_264$zscore,pch=19,cex=0.2,
     xlab="B64_RNAseq_logFC",ylab="TIDE_z_score")
dev.off()


## GSEA analysis on the top and bottom z score genes regarding their expression in 204 and 264
data(examplePathways)
data(exampleRanks)
head(names(examplePathways),4)
head(exampleRanks,4)

high_z_set = high_z_mouse$MGI.symbol
low_z_set = low_z_mouse$MGI.symbol
TIDE_sets = list(high_z_set,low_z_set)
names(TIDE_sets) = c("high_z_genes","low_z_genes")

Rank_204 = expr_diff_204$log2FoldChange
names(Rank_204) = expr_diff_204$symbol
Rank_264 = expr_diff_264$log2FoldChange
names(Rank_264) = expr_diff_264$symbol


GSEA_204 = fgsea(pathways = TIDE_sets, 
                  stats = Rank_204,
                  minSize=15,
                  maxSize=500,
                  nperm=10000)
GSEA_264 = fgsea(pathways = TIDE_sets, 
                 stats = Rank_264,
                 minSize=15,
                 maxSize=500,
                 nperm=10000)

write.csv(GSEA_204[order(pval),1:7],"TIDE_dysfunction/2019-06-13_GSEA_204.csv")
write.csv(GSEA_264[order(pval),1:7],"TIDE_dysfunction/2019-06-13_GSEA_264.csv")

pdf("TIDE_dysfunction/2019-06-13_high_z_204.pdf",width=5,height=3)
plotEnrichment(TIDE_sets[[GSEA_204[order(pval), ]$pathway[1]]],Rank_204)+
  labs(title=GSEA_204[order(pval), ]$pathway[1])
dev.off()

pdf("TIDE_dysfunction/2019-06-13_low_z_204.pdf",width=5,height=3)
plotEnrichment(TIDE_sets[[GSEA_204[order(pval), ]$pathway[2]]],Rank_204)+
  labs(title=GSEA_204[order(pval), ]$pathway[2])
dev.off()

pdf("TIDE_dysfunction/2019-06-13_high_z_264.pdf",width=5,height=3)
plotEnrichment(TIDE_sets[[GSEA_264[order(pval), ]$pathway[1]]],Rank_264)+
  labs(title=GSEA_264[order(pval), ]$pathway[1])
dev.off()

pdf("TIDE_dysfunction/2019-06-13_low_z_264.pdf",width=5,height=3)
plotEnrichment(TIDE_sets[[GSEA_264[order(pval), ]$pathway[2]]],Rank_264)+
  labs(title=GSEA_264[order(pval), ]$pathway[2])
dev.off()
