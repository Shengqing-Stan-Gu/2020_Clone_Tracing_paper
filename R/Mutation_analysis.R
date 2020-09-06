library(maftools)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggpubr)
contrl <- c('P_A','P_B')
treat <- c("s204_A","s204_B","s264_A","s264_B")

## Read in maf files
maf.oncomatrix <- lapply(c(contrl,treat),function(x){
  res.maf <- read.maf(maf=paste(x,".single.varscan.snp.base.filter.maf.dbSNP.maf",sep = ""),
                      useAll = TRUE, isTCGA = FALSE,removeDuplicatedVariants = TRUE)
  oncomatrix <- res.maf@data
  return(oncomatrix)
})
names(maf.oncomatrix) <- c(contrl,treat)

## Compare altered reads frequency in each group
SNP_bulk <- lapply(contrl,function(case){
  cf <- 0.1
  case.oncomatrix <- maf.oncomatrix[[case]] %>%
    dplyr::select(c('Hugo_Symbol','Chromosome','Start_Position', 'End_Position', 'Strand', 'Variant_Classification', 'Variant_Type', 'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2','Tumor_Sample_Barcode','HGVSc', 'HGVSp', 'HGVSp_Short', 't_depth', 't_ref_count', 't_alt_count', 'n_depth', 'n_ref_count', 'n_alt_count')) %>%
    subset(t_alt_count != 0) %>%
    mutate(read.freq = t_alt_count/t_ref_count)
  return(case.oncomatrix)
})
dim(SNP_bulk[[1]])
dim(SNP_bulk[[2]])


SNP_clone <- lapply(treat,function(case){
  cf <- 0.1
  case.oncomatrix <- maf.oncomatrix[[case]] %>%
    dplyr::select(c('Hugo_Symbol','Chromosome','Start_Position', 'End_Position', 'Strand', 'Variant_Classification', 'Variant_Type', 'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2','Tumor_Sample_Barcode','HGVSc', 'HGVSp', 'HGVSp_Short', 't_depth', 't_ref_count', 't_alt_count', 'n_depth', 'n_ref_count', 'n_alt_count')) %>%
    subset(t_alt_count != 0) %>%
    mutate(read.freq = t_alt_count/t_ref_count) %>%
    subset(read.freq > cf)
  return(case.oncomatrix)
})
dim(SNP_clone[[1]])
dim(SNP_clone[[2]])
dim(SNP_clone[[3]])
dim(SNP_clone[[4]])

unique(SNP_clone[[1]]$Hugo_Symbol)
unique(SNP_clone[[2]]$Hugo_Symbol)
unique(SNP_clone[[3]]$Hugo_Symbol)
unique(SNP_clone[[4]]$Hugo_Symbol)

line_204_SNP = intersect(unique(SNP_clone[[1]]$Hugo_Symbol), unique(SNP_clone[[2]]$Hugo_Symbol))
line_264_SNP = intersect(unique(SNP_clone[[3]]$Hugo_Symbol), unique(SNP_clone[[4]]$Hugo_Symbol))

line_204_SNP %in% c(SNP_bulk[[1]]$Hugo_Symbol, SNP_bulk[[2]]$Hugo_Symbol)
line_264_SNP %in% c(SNP_bulk[[1]]$Hugo_Symbol, SNP_bulk[[2]]$Hugo_Symbol)

write.csv(SNP_bulk[[1]], "P_A_mutation.csv", row.names=F)
write.csv(SNP_bulk[[2]], "P_B_mutation.csv", row.names=F)
write.csv(SNP_clone[[1]], "204_A_mutation.csv", row.names=F)
write.csv(SNP_clone[[2]], "204_B_mutation.csv", row.names=F)
write.csv(SNP_clone[[3]], "264_A_mutation.csv", row.names=F)
write.csv(SNP_clone[[4]], "264_B_mutation.csv", row.names=F)





