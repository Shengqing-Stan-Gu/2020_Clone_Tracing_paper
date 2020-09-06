library(ggplot2)
library(ggpubr)
library(cowplot)

metaData = read.csv('../data/ResistClone.Samples.csv')

metaData$File = paste0('../data/TRUST/', metaData$Name, '.bam')

metaData$Name = as.character(metaData$Name)
metaData$Name = sub('-','_',metaData$Name)
metaData$Sample = as.character(metaData$Sample)

metaData$Sample[metaData$Sample == 'Responder_good'] = 'Responder'
metaData$Sample[metaData$Sample == 'Responder_poor'] = 'Nonresponder'

## Default: '20171208', '20171219', '20180608'
metaData = metaData[metaData$Batch %in% c('20171208', '20171219', '20180608'), ]

metaData = metaData[metaData$TimePoint != 'Before', ]
metaData = metaData[metaData$SampleType == 'tumor', ]

samples_of_two_times = c('SG_6L', 'SG_10L', 'SG_18L', 'SG_5R', 'SG_6R', 'SG_7R', 
                         'SG_10R', 'SG_11R', 'SG_14R', 'SG_16R', 'SG_18R', 'SG_19R')
metaData = metaData[!metaData$Name %in% samples_of_two_times, ]



read_trust <- function(f) {
  t = c()
  
  for(suffix in c('-TCR-ALL', '-BCR-HeavyChain', '-BCR-LightChain')) {
    infile = paste0(f, suffix, '.txt')
    
    if(!file.exists(infile)) next ## ignore failed cases
    tab = read.table(infile, sep='\t', header=T)
    if(nrow(tab) == 0) next ## ignore failed cases
    
    tmp = tab[, c('est_lib_size', 'est_clonal_exp', 'cdr3dna', 'cdr3aa', 'Vgene', 'Jgene')]
    colnames(tmp) = c('total', 'count', 'cdr3nt', 'cdr3aa', 'v', 'j')
    tmp$v = sub('([^*]+)[*].+', '\\1', tmp$v)
    tmp$j = sub('([^*]+)[*].+', '\\1', tmp$j)
    tmp$is_complete = grepl('^C[A-Z]*[FW]$', tmp$cdr3aa)
    
    if (suffix == '-BCR-HeavyChain') {
      tmp$c = paste(tab$Cgene, tab$reportgene)
    } else {
      tmp$c = ''
    }
    
    t = rbind(t, tmp)
  }
  
  if(length(t) == 0)
    return(t)
  
  t$v[is.na(t$v)] = ''
  t$j[is.na(t$j)] = ''
  
  t$chain = sapply(t$v, function(x) substr(x,1,3))
  t$chain[t$chain == ''] = sapply(t$j[t$chain == ''], function(x) substr(x,1,3))
  t$chain = unlist(t$chain)
  t$chain[t$chain == 'IGK'] = 'IGL' ## Ig light chain
  
  t = aggregate(count ~ ., t, sum)
  t$freq = 0
  for(c in unique(t$chain))
    t$freq[t$chain == c] = t$count[t$chain == c] / sum(t$count[t$chain == c])
  t = t[order(t$freq, decreasing = T), ]
  
  return(t)
}

data = lapply(as.character(metaData$File), read_trust)
names(data) = metaData$Name


sample_count = as.data.frame.matrix(table(metaData$Batch, metaData$Sample))
print(sample_count)

mouse_site = as.data.frame.matrix(table(metaData$Mouse, metaData$TumorSide))
print(c(sum(mouse_site$Left), sum(mouse_site$Right), sum(mouse_site$Left*mouse_site$Right)))

clone_count = unlist(lapply(data, nrow))
metaData$clone_count = 0
metaData$clone_count[match(names(clone_count), metaData$Name)] = clone_count

ggplot(metaData, aes(Sample, clone_count, color=SampleType)) +
  geom_boxplot() +
  coord_flip() +
  ylab('Number of T and B cell CDR3s') +
  theme_pubr()


get_metric <- function(x, only_complete=T) {
  chains =  c('TRB','TRA','IGH','IGL')
  
  entropy = c()
  clonality = c()
  total = c()
  cpk = c()
  igg = c()
  iga = c()
  for(c in chains) {
    s = x[x$chain == c, ]
    if(length(s) == 0) {
      entropy = c(entropy , NA)
      clonality = c(clonality, NA)
      total = c(total, NA)
      cpk = c(cpk, NA)
      igg = c(igg, NA)
      iga = c(iga, NA)
      next
    }
    if(nrow(s) < 2) {
      entropy = c(entropy , NA)
      clonality = c(clonality, NA)
      total = c(total, max(s$total))
      cpk = c(cpk, length(s) * 1000 / max(s$total))
      igg = c(igg, NA)
      iga = c(iga, NA)
      next
    }
    if(only_complete)
      a = aggregate(freq ~ cdr3nt + cdr3aa + v + j, s[s$is_complete,], sum)
    else
      a = aggregate(freq ~ cdr3nt + cdr3aa + v + j, s, sum)
    p = a$freq
    entropy = c(entropy, sum(-p*log2(p)))
    clonality = c(clonality, 1 - sum(-p*log2(p)) / log2(length(p)))
    total = c(total, max(s$total))
    cpk = c(cpk, length(p) * 1000 / max(s$total))
    
    ig = sum(s$freq) #[grepl('IGH[MDGEA]', s$c)])
    if(ig > 0) {
      igg = c(igg, sum(s$freq[grepl('IGHG', s$c)]) / ig)
      iga = c(iga, sum(s$freq[grepl('IGHA', s$c)]) / ig)
    } else {
      igg = c(igg, NA)
      iga = c(iga, NA)
    }
  }
  return(data.frame(Chain=chains, Entropy=entropy, Clonality=clonality, Total=total, CPK=cpk, IgG=igg, IgA=iga))
}


metrics = do.call(rbind, lapply(data, function(x) get_metric(x, only_complete=TRUE)))
metrics$Name = gsub('[.].*', '', rownames(metrics))
metrics = merge(metaData, metrics, by='Name')
metrics$Chain = factor(metrics$Chain, levels = c('TRB','TRA','IGH','IGL'))

starlog = read.csv('../data/starLog.final.csv')
starlog$Mapped = starlog$Uniquely.mapped.reads.number + starlog$Number.of.reads.mapped.to.multiple.loci +
  starlog$Number.of.reads.mapped.to.too.many.loci
metrics$MappedReads = starlog$Mapped[match(metrics$Name, starlog$Sample)]
metrics$FracRepReads = metrics$Total / metrics$MappedReads

P1 <- ggplot(metrics, aes(Sample, 100*FracRepReads)) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="pointrange") +
  geom_jitter(aes(color=TimePoint, shape=SampleType), width = 0.2) +
  stat_compare_means(aes(label = "p.format"), comparisons = list(c(2,3)), method='wilcox.test') +
  xlab('') +
  ylab('% mapped reads from receptors') +
  scale_y_continuous(limits = c(0, 0.03)) +
  facet_grid(.~Chain) +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none')

P2 <- ggplot(metrics, aes(Sample, Entropy)) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="pointrange") +
  geom_jitter(aes(color=TimePoint, shape=SampleType), width = 0.2) +
  stat_compare_means(aes(label = "p.format"), comparisons = list(c(2,3)), method='wilcox.test') +
  xlab('') +
  ylab("Diversity (Shannon entropy)") +
  scale_y_continuous(limits = c(0, 8)) +
  facet_grid(.~Chain) +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none')

P1
pdf('../work/plot_repertoire_fraction.pdf', height = 4, width = 6)
P1
dev.off()

P2
pdf('../work/plot_repertoire_entropy.pdf', height = 4, width = 6)
P2
dev.off()



