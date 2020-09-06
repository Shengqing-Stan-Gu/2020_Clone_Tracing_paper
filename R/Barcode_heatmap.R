rm(list=ls())

getwd()
library(gplots)

sampleindex = 1:64
samplelist = list()
cumulative = list()

## Read in clontracer output, store as samplelist
for (i in sampleindex) {
  filepath = paste("./2017-09-25_Clontracer_output/New",sprintf("%02d",i),".barcode.txt",sep="")
  samplelist[[i]] = read.table(filepath,header=T,row.names=1,sep="\t",stringsAsFactors=F)
  statpath = paste("./2017-09-25_Clontracer_output/New",sprintf("%02d",i),".cumulative_wealth.txt",sep="")
  cumulative[[i]] = read.table(statpath,header=T,row.names=1,sep="\t",stringsAsFactors=F)
}


###########################
## Analyze in vivo samples
###########################

## Select barcodes in top 90% for each in vivo sample and combine
barcodeid = vector("character")
for (i in 5:64) {
  if(nrow(samplelist[[i]])<10) {
    next
  }
  barcode = rownames(samplelist[[i]])[1:cumulative[[i]][91,1]]
  barcodeid = union(barcode,barcodeid)
}

## Construct data frame with barcodes and frequencies, store as barinvivo.data.frame
barinvivo.data.frame = samplelist[[5]][barcodeid,"Fraction",drop=F]
for (i in 6:length(sampleindex)) {
  barinvivo.data.frame = cbind(barinvivo.data.frame,samplelist[[i]][barcodeid,"Fraction",drop=F])
}
min(barinvivo.data.frame, na.rm=T)
colnames(barinvivo.data.frame) = sprintf("%02d",5:64)
barinvivo.data.frame[c("TCTGAGTCTGAGACAGTGAGAGTGTCTCTC","TGTCAGAGTGACTGAGTGTCTCAGTGAGTC"),]
rownames(barinvivo.data.frame) = barcodeid

## Replace "NA" with 0 and log-transform
barinvivo.data.frame[is.na(barinvivo.data.frame)] = 0
barinvivo.data.frame = barinvivo.data.frame + 1e-6
barinvivo.data.frame = log10(barinvivo.data.frame)
Colsums = colSums(10^(barinvivo.data.frame)-1e-6)

## Remove samples with no read
barinvivo.data.frame = barinvivo.data.frame[,Colsums>0.5]

## Assign response groups for each drug
ctrlsize = c(1258,1066,821,1866,666,859,1066,280,175,
             1680,974,1563,576,2454,1038,1613,1198,2106,951)
pd1size = c(1351,698,543,1474,1103,708,136,86,436,528,97,324,583,336,
            473,374,684,869,287,203,502,507,1400,702,2029,1103,510,260)
ctla4size = c(190,374,1436,1363,187,313,549,116,637,372,359)


#############################################
## Restrict in vivo samples based on response
#############################################

## Select samples with very good or bad response (25 percentiles)
ctrl = which(ctrlsize>=951 & ctrlsize<=1613)
pd1good = which(pd1size<=324)+19
pd1bad = which(pd1size>=869)+19
ctla4good = which(ctla4size<=324)+47
ctla4bad = which(ctla4size>=869)+47

sampleindex = c(ctrl,pd1good,pd1bad,ctla4good,ctla4bad)
barinvivo.data.frame.filter = barinvivo.data.frame[,sampleindex]

## Hierarchical clustering of barinvivo.data.frame.filter
hrow_ward = hclust(dist(barinvivo.data.frame.filter),method="ward")
hcol_ward = hclust(dist(t(barinvivo.data.frame.filter)),method="ward")


## Try making two different heatmaps using different grouping
group1 = c(rep("darkgrey",length(ctrl)),
           rep("red",length(pd1good)+length(pd1bad)),
           rep("orange",length(ctla4good)+length(ctla4bad)))
group2 = c(rep("darkgrey",length(ctrl)),
           rep("blue",length(pd1good)),
           rep("green",length(pd1bad)),
           rep("blue",length(ctla4good)),
           rep("green",length(ctla4bad)))

## Center each barcode to 0
barinvivo.data.frame.filter.norm = barinvivo.data.frame.filter - rowMeans(barinvivo.data.frame.filter)


color.palette = colorRampPalette(c("midnightblue","midnightblue","dodgerblue3","dodgerblue","white","lightpink","indianred1","red2","red2"),space="Lab")
palette.breaks = seq(min(barinvivo.data.frame.filter),max(barinvivo.data.frame.filter),length.out=401)


## Reorder the barinvivo.data.frame.filter.norm by the difference between responder and control group
res_con_diff = rowMeans(barinvivo.data.frame.filter.norm[,c(10:16,24:27)]) - rowMeans(barinvivo.data.frame.filter.norm[,1:7])
barinvivo.data.frame.filter.norm = barinvivo.data.frame.filter.norm[order(res_con_diff),]


pdf("2020-07-21_Heatmap_filter1.pdf",width=14,height=10)
heatmap_filter1 = heatmap.2(as.matrix(barinvivo.data.frame.filter.norm),trace="none",density="none",symm=F,symkey=F,symbreaks=T,
                           Rowv=F,Colv=as.dendrogram(hcol_ward),
                           key=T, scale="row", dendrogram="column", col=color.palette,
                           ColSideColors=group1)
dev.off()

pdf("2020-07-21_Heatmap_filter2.pdf",width=14,height=10)
heatmap_filter2 = heatmap.2(as.matrix(barinvivo.data.frame.filter.norm),trace="none",density="none",symm=F,symkey=F,symbreaks=T,
                            Rowv=F,Colv=as.dendrogram(hcol_ward),
                            key=T, scale="row", dendrogram="column", col=color.palette,
                            ColSideColors=group2)
dev.off()


labrow = rep("", times=dim(barinvivo.data.frame.filter.norm)[1])
labrow[which(rownames(barinvivo.data.frame.filter.norm)=="ACTGAGTGTGAGTGAGTGTCACTGTGAGTC")] = "__________________"
labrow[which(rownames(barinvivo.data.frame.filter.norm)=="TCTGAGTCTGAGACAGTGAGAGTGTCTCTC")] = "__________________"

pdf("2020-07-21_Heatmap_filter1_B04_B64.pdf",width=14,height=20)
heatmap_filter1 = heatmap.2(as.matrix(barinvivo.data.frame.filter.norm),trace="none",density="none",symm=F,symkey=F,symbreaks=T,
                            Rowv=F,Colv=as.dendrogram(hcol_ward),
                            key=T, scale="row", dendrogram="column", col=color.palette,
                            ColSideColors=group1, labRow=labrow)
dev.off()

pdf("2020-07-21_Heatmap_filter1_only_B04_B64.pdf",width=14,height=4)
heatmap_filter1 = heatmap.2(as.matrix(barinvivo.data.frame.filter[c("ACTGAGTGTGAGTGAGTGTCACTGTGAGTC","TCTGAGTCTGAGACAGTGAGAGTGTCTCTC"),]),
                            trace="none",density="none",symm=F,symkey=F,symbreaks=T,
                            Colv=as.dendrogram(hcol_ward),
                            key=T, scale="row", dendrogram="column", col=color.palette,
                            ColSideColors=group1, labRow=c("B04","B64"))
dev.off()

pdf("2020-07-21_Heatmap_filter2_only_B04_B64.pdf",width=14,height=4)
heatmap_filter1 = heatmap.2(as.matrix(barinvivo.data.frame.filter[c("ACTGAGTGTGAGTGAGTGTCACTGTGAGTC","TCTGAGTCTGAGACAGTGAGAGTGTCTCTC"),]),
                            trace="none",density="none",symm=F,symkey=F,symbreaks=T,
                            Colv=as.dendrogram(hcol_ward),
                            key=T, scale="row", dendrogram="column", col=color.palette,
                            ColSideColors=group2, labRow=c("B04","B64"))
dev.off()

