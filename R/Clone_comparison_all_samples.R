rm(list=ls())

getwd()

barinvivo.data.frame = read.csv("2017-09-26_Barcode_freq_invivo.csv", header=T, row.names=1, stringsAsFactors=F)

## Select barcodes with mean log10 frequency > -4
barinvivo.data.frame = barinvivo.data.frame[rowMeans(barinvivo.data.frame)>(-4),]


## Define groups
ctrlsize = c(1258,1066,821,1866,666,859,1066,280,175,
             1680,974,1563,576,2454,1038,1613,1198,2106,951)
pd1size = c(1351,698,543,1474,1103,708,136,86,436,528,97,324,583,336,
            473,374,684,869,287,203,502,507,1400,702,2029,1103,510,260)
ctla4size = c(190,374,1436,1363,187,313,549,116,637,372,359)
pd1 = which(pd1size>1)+19
ctla4 = which(ctla4size>1)+47


## Perform t.test between IgG control and anti-PD1 recipients
pd1_ttestpvalues = vector()
pd1_tvalues = vector()
for (i in 1:dim(barinvivo.data.frame)[1]) {
  if (var(as.numeric(barinvivo.data.frame[i,c(1:19,pd1)]))==0) {
    pd1_ttestpvalues = c(pd1_ttestpvalues,1)
    pd1_tvalues = c(pd1_tvalues,0)
  }
  else {
    pd1_ttestresult = t.test(barinvivo.data.frame[i,1:19], barinvivo.data.frame[i,pd1])
    pd1_ttestpvalues = c(pd1_ttestpvalues, pd1_ttestresult$p.value)
    pd1_tvalues = c(pd1_tvalues, pd1_ttestresult$statistic)
  }
}

pd1_pvaluerank = rank(pd1_ttestpvalues)
pd1_ttestFDR = pd1_ttestpvalues*length(pd1_ttestpvalues)/pd1_pvaluerank
write.csv(cbind(barinvivo.data.frame,pd1_ttestpvalues,pd1_ttestFDR,pd1_tvalues)[order(pd1_ttestFDR),],
          "./2018-11-20_IgG_PD1_comparison.csv")



## Perform t.test between IgG control and anti-CTLA4 recipients
ctla4_ttestpvalues = vector()
ctla4_tvalues = vector()
for (i in 1:dim(barinvivo.data.frame)[1]) {
  if (var(as.numeric(barinvivo.data.frame[i,c(1:19,ctla4)]))==0) {
    ctla4_ttestpvalues = c(ctla4_ttestpvalues,1)
    ctla4_tvalues = c(ctla4_tvalues,0)
  }
  else {
    ctla4_ttestresult = t.test(barinvivo.data.frame[i,1:19], barinvivo.data.frame[i,ctla4])
    ctla4_ttestpvalues = c(ctla4_ttestpvalues,ctla4_ttestresult$p.value)
    ctla4_tvalues = c(ctla4_tvalues,ctla4_ttestresult$statistic)
  }
}

ctla4_pvaluerank = rank(ctla4_ttestpvalues)
ctla4_ttestFDR = ctla4_ttestpvalues*length(ctla4_ttestpvalues)/ctla4_pvaluerank
write.csv(cbind(barinvivo.data.frame,ctla4_ttestpvalues,ctla4_ttestFDR,ctla4_tvalues)[order(ctla4_ttestFDR),],
          "./2018-11-20_IgG_CTLA4_comparison.csv")


## Compare the enrichment under anti-PD1 and anti-CTLA4
pd1_ctla4_cor = cor.test(pd1_tvalues,ctla4_tvalues,method="spearman")
pd1_ctla4_cor

library(ggplot2)
Cor_plot = qplot(x=-pd1_tvalues,y=-ctla4_tvalues)+
  geom_point(color="black",size=0.8)+
  xlab("t values of anti-PD1 vs control")+
  ylab("t values of anti-CTLA4 vs control")+
  theme_bw()+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14))
pdf("./2018-11-26_PD1_CTLA4_correlation.pdf",width=5,height=4)
Cor_plot
dev.off()


## Compare the overlap between anti-PD1 and anti-CTLA4
pd1_enriched = rownames(barinvivo.data.frame)[pd1_ttestFDR<0.2 & pd1_tvalues<0]
pd1_depleted = rownames(barinvivo.data.frame)[pd1_ttestFDR<0.2 & pd1_tvalues>0]
ctla4_enriched = rownames(barinvivo.data.frame)[ctla4_ttestFDR<0.2 & ctla4_tvalues<0]
ctla4_depleted = rownames(barinvivo.data.frame)[ctla4_ttestFDR<0.2 & ctla4_tvalues>0]
intersect(pd1_enriched,ctla4_enriched)
intersect(pd1_depleted,ctla4_depleted)
enriched = union(pd1_enriched, ctla4_enriched)
depleted = union(pd1_depleted, ctla4_depleted)


## Compare the percentages of enriched clones in different groups
enriched.data.frame = barinvivo.data.frame[enriched,]
depleted.data.frame = barinvivo.data.frame[depleted,]
enriched.data.frame[1:4,1:4]
enriched.data.frame = 10^enriched.data.frame - 10^(-6)
depleted.data.frame = 10^depleted.data.frame - 10^(-6)
enriched_sum_freq = colSums(enriched.data.frame)
depleted_sum_freq = colSums(depleted.data.frame)
enriched_sum_freq
identity = c(rep("IgG",19),rep("PD1",length(pd1)),rep("CTLA4",length(ctla4)))
write.csv(cbind(identity,enriched_sum_freq), "2018-11-26_Resistant_frequencies.csv")
write.csv(cbind(identity,depleted_sum_freq), "2018-11-26_Sensitive_frequencies.csv")


## Compare the percentages of line 204 and 264 in different groups
freq.204 = 10^barinvivo.data.frame["ACTGAGTGTGAGTGAGTGTCACTGTGAGTC",] - 10^(-6)
freq.264 = 10^barinvivo.data.frame["TCTGAGTCTGAGACAGTGAGAGTGTCTCTC",] - 10^(-6)
write.csv(cbind(identity,t(freq.204),t(freq.264)), "2018-11-27_Select_line_frequencies.csv")



