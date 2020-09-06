library("EnsDb.Mmusculus.v79")
library("EnsDb.Hsapiens.v79")
library("org.Hs.eg.db")
library("annotate")
library("biomaRt")
library("ggplot2")
library(survival)
library(survminer)

getwd()

##########################################
### Normalize clinical expression data ###
##########################################

Expr_pre = read.table("../2019-06-29_Line_204_264_in_clinical/ICB_Clinical_data/Mariathasan2018_PDL1/IMvigor210.self_subtract",sep="\t",row.names=1,stringsAsFactors=F)
Expr_pre[1:4,1:4]

Outcome = read.table("../2019-06-29_Line_204_264_in_clinical/ICB_Clinical_data/Mariathasan2018_PDL1/IMvigor210.clinical",sep="\t",header=T,row.names=1,stringsAsFactors=F)
head(Outcome,4)

table(colnames(Expr_pre)==rownames(Outcome))

summary(as.factor(Outcome$Tissue))
bladder_index = which(Outcome$Tissue=="bladder")
Expr_pre = Expr_pre[,bladder_index]
Outcome = Outcome[bladder_index,]

table(colnames(Expr_pre)==rownames(Outcome))

##################################################
### Convert human gene ID to mouse gene symbol ###
##################################################

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
Pre_symbol_mouse = getLDS(attributes = "entrezgene_id",
                          filters = "entrezgene_id",
                          values = rownames(Expr_pre),
                          mart = human,
                          attributesL = "mgi_symbol",
                          martL = mouse,
                          uniqueRows = T)
head(Pre_symbol_mouse,4)
Pre_symbol_mouse = Pre_symbol_mouse[!duplicated(Pre_symbol_mouse$MGI.symbol),]
Pre_symbol_mouse = Pre_symbol_mouse[!duplicated(Pre_symbol_mouse$NCBI.gene.ID),]


shared_index = vector()
shared_symbol = data.frame()
for (i in 1:dim(Expr_pre)[1]) {
  if (rownames(Expr_pre[i,]) %in% Pre_symbol_mouse$NCBI.gene.ID) {
    shared_index = c(shared_index, i)
    shared_symbol =  rbind(shared_symbol, Pre_symbol_mouse[which(Pre_symbol_mouse$NCBI.gene.ID==rownames(Expr_pre[i,])),])
  }
}
Expr_mgi = Expr_pre[shared_index,]
table(rownames(Expr_mgi)==shared_symbol$NCBI.gene.ID)
rownames(Expr_mgi) = shared_symbol$MGI.symbol


###########################################
### Get mouse gene expression signature ###
###########################################

Deseq_out = read.csv("2018-12-05_Differential_Expression_poor_good.csv",row.names=1,stringsAsFactors=F)
Deseq_out = Deseq_out[which(Deseq_out$padj<0.1),]
Deseq_out = Deseq_out[order(abs(Deseq_out$log2FoldChange), decreasing=T),]
head(Deseq_out,4)
## Convert log fold change to responder minus nonresponder
Deseq_out$log2FoldChange = -Deseq_out$log2FoldChange
head(Deseq_out,4)

Sig_invivo = Deseq_out[1:200, c("symbol","log2FoldChange")]
Sig_invivo = Sig_invivo[which(Sig_invivo$symbol %in% rownames(Expr_mgi)),]


####################################################################
### Compare responders and non-resonders in resistance signature ###
####################################################################

Tumor_invivo = Expr_mgi[Sig_invivo$symbol,]

cor_invivo = vector("numeric")
for (i in 1:dim(Tumor_invivo)[2]) {
  cor_invivo = c(cor_invivo, cor(Tumor_invivo[,i], Sig_invivo[,2]))
}
names(cor_invivo) = colnames(Tumor_invivo)


## Separate correlations between responders and non-responders and compare correlations
Res_all = which(Outcome$Best %in% c(3))
Nonres_all = which(Outcome$Best %in% c(0))

t.test(cor_invivo[Res_all], cor_invivo[Nonres_all])

write.csv(cor_invivo[Res_all], "Exploratory_analysis//2019-07-17_Mariathasan_CR_pearson.csv")
write.csv(cor_invivo[Nonres_all], "Exploratory_analysis/2019-07-17_Mariathasan_PD_pearson.csv")


#######################################################################
### Draw survival curve based on correlation with in vivo signature ###
#######################################################################

summary(cor_invivo)
cor_order = order(cor_invivo, decreasing=T)
patient_survival = cbind(cor_invivo[cor_order],Outcome[cor_order,c(1:2,4)])
write.csv(patient_survival, "Exploratory_analysis/2019-07-17_Mariathasan_survival.csv",row.names=F)





