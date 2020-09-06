library(ggplot2)
library(flexmix)

set.seed(100)

barcode_freq = read.csv('2019-02-27_Barcode_freq_invivo_freq.csv', stringsAsFactors = F, row.names = 1)
barcode_sample_info = read.csv(("2019-02-27_Barcode_freq_invivo_sample_info.csv"), stringsAsFactors=F, row.names = 1)

index_present = which(rowMeans(barcode_freq[,1:19])>(-4))
barcode_freq = barcode_freq[index_present,]

barcode_freq_control = rowMeans(barcode_freq[,1:19])
barcode_freq_PD1 = rowMeans(barcode_freq[,20:47])
barcode_freq_CTLA4 = rowMeans(barcode_freq[,48:58])
group_means = cbind(barcode_freq_control, barcode_freq_PD1, barcode_freq_CTLA4)

barcode_means = rowMeans(barcode_freq)

bfreq = NULL
N = nrow(barcode_freq)
for(i in colnames(barcode_freq)) {
  val = barcode_freq[[i]]
  info = barcode_sample_info[[i]]
  new = data.frame(Sample = i, 
                   Barcode = rownames(barcode_freq), 
                   Treatment = info[1],
                   Tumor_size = as.numeric(info[2]),
                   Response_assignment = info[3],
                   RNAseq_available = info[4],
                   Freq = as.numeric(val[1:N]),
                   Raw = 10**as.numeric(val[1:N]) - 1e-6)
  bfreq = rbind(bfreq, new)
}



##ggplot(bfreq, aes(x=Freq)) + 
##  geom_histogram(bins = 50, aes(y=..density..), position="identity", alpha=0.5) +
##  geom_density(alpha=0.6) +
##  facet_wrap(~Sample+Response_assignment, nrow=7, scales="free_y")
##dev.off()


bfreq$Treatment = factor(as.character(bfreq$Treatment), levels=c('control_IgG', 'anti-PD1', 'anti-CTLA4'))
bfreq$Freq = bfreq$Freq - barcode_freq_control

mm = stepFlexmix(Freq ~ Treatment | Barcode, data = bfreq, 
                 k = 1:6, nrep = 3)

m2 = getModel(mm, "BIC")
summary(m2)

rm1 = refit(m2)
summary(rm1)


pdf("2019-06-18_Cluster_assignment.pdf", height = 4, width = 8)
plot(m2)
dev.off()

psize = prior(m2)
coefs = parameters(m2)

responseorder = order(coefs[3,])
psize = psize[responseorder]
coefs = coefs[,responseorder]


mpar = NULL
for(i in 1:length(psize)) {
  model_par = data.frame(Cluster=i, Parameter=sub('Treatment','',row.names(coefs)), Value=coefs[,i])
  model_par = rbind(data.frame(Cluster=i, Parameter='ClusterSize', Value=psize[i]), model_par)
  mpar = rbind(mpar, model_par)
}
mpar = mpar[mpar$Parameter != 'sigma', ]
mpar = mpar[mpar$Parameter != 'coef.(Intercept)',]
mpar$Cluster = factor(mpar$Cluster)
##levels(mpar$Cluster) = length(psize) - order(psize) + 1 ## order clusters by size
mpar$Cluster = as.character(mpar$Cluster)

mpar = mpar[row.names(mpar)!="coef.(Intercept)*",]

pdf("2019-06-18_Flexmix_summary.pdf", height = 3, width = 8)
ggplot(mpar, aes(Cluster, Value)) +
  geom_bar(stat="identity") +
  facet_wrap(~Parameter, nrow=1, scales="free_y") +
  theme_bw()
dev.off()

pdf("2019-06-18_Cluster_size.pdf", height=3, width=3)
ggplot(mpar[mpar$Parameter=="ClusterSize",], aes(Cluster, Value)) +
  geom_bar(stat="identity") +
  facet_wrap(~Parameter, scales="free_y") +
  theme_bw()
dev.off()

pdf("2019-06-18_PD1_coef.pdf", height=3, width=3)
ggplot(mpar[mpar$Parameter=="coef.anti-PD1",], aes(Cluster, Value)) +
  geom_bar(stat="identity", fill="darkgreen") +
  facet_wrap(~Parameter, scales="free_y") +
  theme_bw()
dev.off()

pdf("2019-06-18_CTLA4_coef.pdf", height=3, width=3)
ggplot(mpar[mpar$Parameter=="coef.anti-CTLA4",], aes(Cluster, Value)) +
  geom_bar(stat="identity", fill="darkblue") +
  facet_wrap(~Parameter, scales="free_y") +
  theme_bw()
dev.off()


head(posterior(m2),20)

cluster_assignment = posterior(m2)[,responseorder]
for(i in 1:length(psize)) {
  bfreq[[paste0('Cluster',i)]] = cluster_assignment[, i]
}
write.csv(bfreq, file='2019-06-18_barcode_model_fitted.csv', row.names = F)





barcode_assignment = cluster_assignment[1:N,]

cluster1 = which(barcode_assignment[,1]>0.5)
cluster2 = which(barcode_assignment[,2]>0.5)
cluster3 = which(barcode_assignment[,3]>0.5)
cluster4 = which(barcode_assignment[,4]>0.5)
cluster5 = which(barcode_assignment[,5]>0.5)
cluster6 = which(barcode_assignment[,6]>0.5)


barcode_freq_linear = barcode_freq
for (i in 1:dim(barcode_freq)[2]) {
  barcode_freq_linear[,i] = 10^as.numeric(barcode_freq[,i]) - 1E-6
}

Sum_cluster1 = colSums(barcode_freq_linear[cluster1,])
Sum_cluster2 = colSums(barcode_freq_linear[cluster2,])
Sum_cluster3 = colSums(barcode_freq_linear[cluster3,])
Sum_cluster4 = colSums(barcode_freq_linear[cluster4,])
Sum_cluster5 = colSums(barcode_freq_linear[cluster5,])
Sum_cluster6 = colSums(barcode_freq_linear[cluster6,])


Cluster_freq = cbind(Sum_cluster1,Sum_cluster2,Sum_cluster3,Sum_cluster4,Sum_cluster5,Sum_cluster6,
                     size=log10(as.numeric(barcode_sample_info["Tumor_size",])))
write.csv(Cluster_freq, "2019-08-18_Cluster_freq.csv")


Resistant_cluster = Sum_cluster5+Sum_cluster6
plot(Resistant_cluster[20:47], log10(as.numeric(barcode_sample_info["Tumor_size",20:47])))
plot(Resistant_cluster[48:58], log10(as.numeric(barcode_sample_info["Tumor_size",48:58])))
plot(Resistant_cluster[1:19], barcode_sample_info["Tumor_size",1:19])
cor_PD1 = cor.test(Resistant_cluster[20:47], as.numeric(barcode_sample_info["Tumor_size",20:47]), method="spearman")
cor_CTLA4 = cor.test(Resistant_cluster[48:58], as.numeric(barcode_sample_info["Tumor_size",48:58]), method="spearman")
cor_control = cor.test(Resistant_cluster[1:19], as.numeric(barcode_sample_info["Tumor_size",1:19]), method="spearman")
cor_PD1
cor_CTLA4
cor_control
cor.test(Resistant_cluster, as.numeric(barcode_sample_info["Tumor_size",]), method="spearman")


pdf("2019-06-18_Resistant_frequency_control.pdf", height=3, width=4.5)
ggplot(data.frame(Frequency=Resistant_cluster[1:19], Size=log10(as.numeric(barcode_sample_info["Tumor_size",1:19]))),
       aes(x=Frequency, y=Size)) +
  geom_point(color="black") +
  xlab("Frequency") +
  ylab(expression("log"[10]*"(Tumor_size)")) +
  theme_bw()
dev.off()

pdf("2019-06-18_Resistant_frequency_PD1.pdf", height=3, width=4.5)
ggplot(data.frame(Frequency=Resistant_cluster[20:47], Size=log10(as.numeric(barcode_sample_info["Tumor_size",20:47]))),
       aes(x=Frequency, y=Size)) +
  geom_point(color="darkgreen") +
  xlab("Frequency") +
  ylab(expression("log"[10]*"(Tumor_size)")) +
  theme_bw()
dev.off()

pdf("2019-06-18_Resistant_frequency_CTLA4.pdf", height=3, width=4.5)
ggplot(data.frame(Frequency=Resistant_cluster[48:58], Size=log10(as.numeric(barcode_sample_info["Tumor_size",48:58]))),
       aes(x=Frequency, y=Size)) +
  geom_point(color="darkblue") +
  xlab("Frequency") +
  ylab(expression("log"[10]*"(Tumor_size)")) +
  theme_bw()
dev.off()


Sensitive_cluster = Sum_cluster1+Sum_cluster2+Sum_cluster3
cor_PD1 = cor.test(Sensitive_cluster[20:47], as.numeric(barcode_sample_info["Tumor_size",20:47]), method="spearman")
cor_CTLA4 = cor.test(Sensitive_cluster[48:58], as.numeric(barcode_sample_info["Tumor_size",48:58]), method="spearman")
cor_control = cor.test(Sensitive_cluster[1:19], as.numeric(barcode_sample_info["Tumor_size",1:19]), method="spearman")
cor_PD1
cor_CTLA4
cor_control
cor.test(Sensitive_cluster, as.numeric(barcode_sample_info["Tumor_size",]), method="spearman")

pdf("2019-06-18_Sensitive_frequency_control.pdf", height=3, width=4.5)
ggplot(data.frame(Frequency=Sensitive_cluster[1:19], Size=log10(as.numeric(barcode_sample_info["Tumor_size",1:19]))),
       aes(x=Frequency, y=Size)) +
  geom_point(color="black") +
  xlab("Frequency") +
  ylab(expression("log"[10]*"(Tumor_size)")) +
  theme_bw()
dev.off()

pdf("2019-06-18_Sensitive_frequency_PD1.pdf", height=3, width=4.5)
ggplot(data.frame(Frequency=Sensitive_cluster[20:47], Size=log10(as.numeric(barcode_sample_info["Tumor_size",20:47]))),
       aes(x=Frequency, y=Size)) +
  geom_point(color="darkgreen") +
  xlab("Frequency") +
  ylab(expression("log"[10]*"(Tumor_size)")) +
  theme_bw()
dev.off()

pdf("2019-06-18_Sensitive_frequency_CTLA4.pdf", height=3, width=4.5)
ggplot(data.frame(Frequency=Sensitive_cluster[48:58], Size=log10(as.numeric(barcode_sample_info["Tumor_size",48:58]))),
       aes(x=Frequency, y=Size)) +
  geom_point(color="darkblue") +
  xlab("Frequency") +
  ylab(expression("log"[10]*"(Tumor_size)")) +
  theme_bw()
dev.off()

