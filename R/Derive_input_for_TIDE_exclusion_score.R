library(biomaRt)

rm(list=ls())

## Read in the rsem output (TPM)
rsem_tpm = read.csv("./RSEM_result/2019-01-02_Scaled_TPM_with_Symbol.csv",row.names=1,stringsAsFactors=F)
colnames(rsem_tpm) = c("14L","14R","16R","17R","20R","24L","26R","29L",
                         "30L","30R","32R","33L","33R","35L","3_invitro","4_invitro","symbol")
rsem_tpm[1:5,1:5]

response = c("control IgG","control IgG","control IgG","control IgG","control IgG",
             "responder","responder","non-responder","responder","responder","non-responder","non-responder","non-responder","responder")


## Get human gene symbol
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
mouse_to_human = getLDS(attributes = "mgi_symbol",
                        filters = "mgi_symbol",
                        values = rsem_tpm$symbol,
                        mart = mouse,
                        attributesL = "hgnc_symbol",
                        martL = human,
                        uniqueRows = T)
head(mouse_to_human,5)


## Select only the rows of genes with MGI symbol and focus on relative expression
rsem_tpm = rsem_tpm[which(rsem_tpm$symbol %in% mouse_to_human$MGI.symbol),]
relative_expr = rsem_tpm
relative_expr_merged = merge(relative_expr,mouse_to_human,by.x="symbol",by.y="MGI.symbol",sort=F)
relative_expr_merged$symbol = relative_expr_merged$HGNC.symbol
relative_expr_merged = relative_expr_merged[,1:15]
head(relative_expr_merged,4)

relative_expr_merged[,2:15] = log2(relative_expr_merged[,2:15]+1)
mean_expr = rowMeans(relative_expr_merged[,2:15])
relative_expr_merged[,2:15] = relative_expr_merged[,2:15]-mean_expr
head(relative_expr_merged,4)


## Write output
write.table(relative_expr_merged[,2:15],"2019-02-08_Relative_expr",sep="\t",
            row.names=relative_expr_merged[,1],quote=F)
