library("EnsDb.Mmusculus.v79")
library("EnsDb.Hsapiens.v79")
library("biomaRt")

getwd()

Expr_pre = read.table("Riaz/Mel_PD1_Riaz.expression.Pre",sep="\t",row.names=1,stringsAsFactors=F)
colnames(Expr_pre)

Outcome = read.csv("Riaz/2019-07-12_Cleaned_outcome.csv",header=T,row.names=1,stringsAsFactors=F)
head(Outcome,4)

## Calculate the expression of each gene in each sample normalized to the row mean
expr_mean = rowMeans(Expr_pre)
Expr_pre_norm = Expr_pre - expr_mean
head(rowMeans(Expr_pre_norm),4)

## Import the expression signatures of line 204 and 264
Deseq_204 = read.csv("Line204_vs_parental.csv",row.names=1,stringsAsFactors=F)
Deseq_264 = read.csv("Line264_vs_parental.csv",row.names=1,stringsAsFactors=F)

Deseq_204 = Deseq_204[which(Deseq_204$padj<0.1),]
Deseq_264 = Deseq_264[which(Deseq_264$padj<0.1),]

Deseq_204 = Deseq_204[order(abs(Deseq_204$log2FoldChange), decreasing=T),]
Deseq_264 = Deseq_264[order(abs(Deseq_264$log2FoldChange), decreasing=T),]

## Get the signature of line 204 and line 264
Sig_204 = Deseq_204[1:200,c("symbol","log2FoldChange")]
Sig_264 = Deseq_264[1:200,c("symbol","log2FoldChange")]


## Convert human gene symbol in clinical cohort into mouse gene symbol
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
Pre_symbol_mouse = getLDS(attributes = "hgnc_symbol",
                          filters = "hgnc_symbol",
                          values = row.names(Expr_pre_norm),
                          mart = human,
                          attributesL = "mgi_symbol",
                          martL = mouse,
                          uniqueRows = T)
head(Pre_symbol_mouse,4)


## Construct two new data.frames for clinical data and mouse cell line
Expr_pre_204_ordered = data.frame()
Expr_pre_264_ordered = data.frame()
Sig_204_ordered = data.frame()
Sig_264_ordered = data.frame()

Human_genes = row.names(Expr_pre_norm)
Mouse_genes_204 = Sig_204$symbol
Mouse_genes_264 = Sig_264$symbol

for (i in 1:dim(Pre_symbol_mouse)[1]) {
  if (Pre_symbol_mouse$HGNC.symbol[i] %in% Human_genes) {
    if (Pre_symbol_mouse$MGI.symbol[i] %in% Mouse_genes_204) {
      temp_human = Expr_pre_norm[which(Human_genes == Pre_symbol_mouse$HGNC.symbol[i])[1],]
      Expr_pre_204_ordered = rbind(Expr_pre_204_ordered,temp_human)
      
      temp_mouse = Sig_204[which(Sig_204$symbol==Pre_symbol_mouse$MGI.symbol[i])[1],]
      Sig_204_ordered = rbind(Sig_204_ordered,temp_mouse)
    }
  }
}

for (i in 1:dim(Pre_symbol_mouse)[1]) {
  if (Pre_symbol_mouse$HGNC.symbol[i] %in% Human_genes) {
    if (Pre_symbol_mouse$MGI.symbol[i] %in% Mouse_genes_264) {
      temp_human = Expr_pre_norm[which(Human_genes == Pre_symbol_mouse$HGNC.symbol[i])[1],]
      Expr_pre_264_ordered = rbind(Expr_pre_264_ordered,temp_human)
      
      temp_mouse = Sig_264[which(Sig_264$symbol==Pre_symbol_mouse$MGI.symbol[i])[1],]
      Sig_264_ordered = rbind(Sig_264_ordered,temp_mouse)
    }
  }
}

head(Expr_pre_204_ordered, 4)
head(Sig_204_ordered, 4)

## Calculate the pearson correlation between each patient and each resistance signature
Pearson_204 = vector("numeric")
for (i in 1:dim(Expr_pre_204_ordered)[2]) {
  Pearson_204 = c(Pearson_204, cor(Expr_pre_204_ordered[,i], Sig_204_ordered[,2], method="pearson"))
}
names(Pearson_204) = colnames(Expr_pre_204_ordered)
Pearson_204

Pearson_264 = vector("numeric")
for (i in 1:dim(Expr_pre_264_ordered)[2]) {
  Pearson_264 = c(Pearson_264, cor(Expr_pre_264_ordered[,i], Sig_264_ordered[,2], method="pearson"))
}
names(Pearson_264) = colnames(Expr_pre_264_ordered)
Pearson_264

table(names(Pearson_204) == names(Pearson_264))


## Separate correlations between responders and non-responders and compare correlations
Res_204 = vector()
Res_264 = vector()
Nonres_204 = vector()
Nonres_264 = vector()

for (i in 1:length(Pearson_204)) {
  if (Outcome[names(Pearson_204[i]),"Outcome"] %in% c("PRCR")) {
    Res_204 = c(Res_204,Pearson_204[i])
    Res_264 = c(Res_264,Pearson_264[i])
  }
  else if (Outcome[names(Pearson_204[i]),"Outcome"] %in% c("PD")) {
    Nonres_204 = c(Nonres_204,Pearson_204[i])
    Nonres_264 = c(Nonres_264,Pearson_264[i])
  }
}

t.test(Res_204, Nonres_204)
t.test(Res_264, Nonres_264)



write.csv(cbind(Res_204,Res_264), "2019-07-11_Responder_pretreatment_pearson.csv")
write.csv(cbind(Nonres_204,Nonres_264), "2019-07-11_Nonresponder_pretreatment_pearson.csv")






Expr_pre = read.table("Riaz/Mel_PD1_Riaz.expression.On",sep="\t",row.names=1,stringsAsFactors=F)
colnames(Expr_pre)

Outcome = read.csv("Riaz/2019-07-12_Cleaned_outcome.csv",header=T,row.names=1,stringsAsFactors=F)
head(Outcome,4)

## Calculate the expression of each gene in each sample normalized to the row mean
expr_mean = rowMeans(Expr_pre)
Expr_pre_norm = Expr_pre - expr_mean
head(rowMeans(Expr_pre_norm),4)

## Import the expression signatures of line 204 and 264
Deseq_204 = read.csv("Line204_vs_parental.csv",row.names=1,stringsAsFactors=F)
Deseq_264 = read.csv("Line264_vs_parental.csv",row.names=1,stringsAsFactors=F)

Deseq_204 = Deseq_204[which(Deseq_204$padj<0.1),]
Deseq_264 = Deseq_264[which(Deseq_264$padj<0.1),]

Deseq_204 = Deseq_204[order(abs(Deseq_204$log2FoldChange), decreasing=T),]
Deseq_264 = Deseq_264[order(abs(Deseq_264$log2FoldChange), decreasing=T),]

## Get the signature of line 204 and line 264
Sig_204 = Deseq_204[1:200,c("symbol","log2FoldChange")]
Sig_264 = Deseq_264[1:200,c("symbol","log2FoldChange")]


## Convert human gene symbol in clinical cohort into mouse gene symbol
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
Pre_symbol_mouse = getLDS(attributes = "hgnc_symbol",
                      filters = "hgnc_symbol",
                      values = row.names(Expr_pre_norm),
                      mart = human,
                      attributesL = "mgi_symbol",
                      martL = mouse,
                      uniqueRows = T)
head(Pre_symbol_mouse,4)


## Construct two new data.frames for clinical data and mouse cell line
Expr_pre_204_ordered = data.frame()
Expr_pre_264_ordered = data.frame()
Sig_204_ordered = data.frame()
Sig_264_ordered = data.frame()

Human_genes = row.names(Expr_pre_norm)
Mouse_genes_204 = Sig_204$symbol
Mouse_genes_264 = Sig_264$symbol

for (i in 1:dim(Pre_symbol_mouse)[1]) {
  if (Pre_symbol_mouse$HGNC.symbol[i] %in% Human_genes) {
    if (Pre_symbol_mouse$MGI.symbol[i] %in% Mouse_genes_204) {
      temp_human = Expr_pre_norm[which(Human_genes == Pre_symbol_mouse$HGNC.symbol[i])[1],]
      Expr_pre_204_ordered = rbind(Expr_pre_204_ordered,temp_human)
      
      temp_mouse = Sig_204[which(Sig_204$symbol==Pre_symbol_mouse$MGI.symbol[i])[1],]
      Sig_204_ordered = rbind(Sig_204_ordered,temp_mouse)
    }
  }
}

for (i in 1:dim(Pre_symbol_mouse)[1]) {
  if (Pre_symbol_mouse$HGNC.symbol[i] %in% Human_genes) {
    if (Pre_symbol_mouse$MGI.symbol[i] %in% Mouse_genes_264) {
      temp_human = Expr_pre_norm[which(Human_genes == Pre_symbol_mouse$HGNC.symbol[i])[1],]
      Expr_pre_264_ordered = rbind(Expr_pre_264_ordered,temp_human)
      
      temp_mouse = Sig_264[which(Sig_264$symbol==Pre_symbol_mouse$MGI.symbol[i])[1],]
      Sig_264_ordered = rbind(Sig_264_ordered,temp_mouse)
    }
  }
}

head(Expr_pre_204_ordered, 4)
head(Sig_204_ordered, 4)
head(Expr_pre_264_ordered, 4)
head(Sig_264_ordered, 4)

## Calculate the pearson correlation between each patient and each resistance signature
Pearson_204 = vector("numeric")
for (i in 1:dim(Expr_pre_204_ordered)[2]) {
  Pearson_204 = c(Pearson_204, cor(Expr_pre_204_ordered[,i], Sig_204_ordered[,2], method="pearson"))
}
names(Pearson_204) = colnames(Expr_pre_204_ordered)
Pearson_204

Pearson_264 = vector("numeric")
for (i in 1:dim(Expr_pre_264_ordered)[2]) {
  Pearson_264 = c(Pearson_264, cor(Expr_pre_264_ordered[,i], Sig_264_ordered[,2], method="pearson"))
}
names(Pearson_264) = colnames(Expr_pre_264_ordered)
Pearson_264

table(names(Pearson_204) == names(Pearson_264))


## Separate correlations between responders and non-responders and compare correlations
Res_204 = vector()
Res_264 = vector()
Nonres_204 = vector()
Nonres_264 = vector()

for (i in 1:length(Pearson_204)) {
  if (Outcome[names(Pearson_204[i]),"Outcome"] %in% c("PRCR")) {
    Res_204 = c(Res_204,Pearson_204[i])
    Res_264 = c(Res_264,Pearson_264[i])
  }
  else if (Outcome[names(Pearson_204[i]),"Outcome"] %in% c("PD")) {
    Nonres_204 = c(Nonres_204,Pearson_204[i])
    Nonres_264 = c(Nonres_264,Pearson_264[i])
  }
}

t.test(Res_204, Nonres_204)
t.test(Res_264, Nonres_264)



write.csv(cbind(Res_204,Res_264), "2019-07-11_Responder_ontreatment_pearson.csv")
write.csv(cbind(Nonres_204,Nonres_264), "2019-07-11_Nonresponder_ontreatment_pearson.csv")


