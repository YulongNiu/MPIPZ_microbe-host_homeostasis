######################hierarchical clustering####################
setwd('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/')

library('readr')
library('magrittr')
library('tibble')
library('gplots')
library('dendextend')
library('dynamicTreeCut')
library('ggplot2')
library('tidyr')
library('DESeq2')
library('dplyr')
library('RColorBrewer')

load('degres_condi_Mock.RData')
degres <- read_csv('eachGroup_vs_Mock_k.csv',
                   col_types = cols(Chromosome = col_character()))

##~~~~~~~~~~~~~~~~~~~~~~useful funcs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
meanFlg22 <- function(v) {

  require('magrittr')

  res <- v %>%
    split(rep(1 : 4, each = 3)) %>%
    sapply(mean, na.rm = TRUE)

  return(res)
}

##p value calculation from WGCNA
corPvalueStudent <- function(cor, nSamples) {

  ## ref: https://courses.lumenlearning.com/introstats1/chapter/testing-the-significance-of-the-correlation-coefficient/
  T <- sqrt(nSamples - 2) * cor / sqrt(1 - cor^2)

  p <- apply(T, 1:2, function(x) {
    if (x < 0) {
      eachp <- 1 -  pt(x, nSamples - 2, lower.tail = FALSE)
    } else {
      eachp <- pt(x, nSamples - 2, lower.tail = FALSE)
    }
  })

  return(p)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~prepare counts~~~~~~~~~~~~~~~~~~~~~~~~~~
rawCount <- counts(degres)

## mean value of normalized count
meanCount <- rawCount %>%
  apply(1, meanFlg22) %>%
  t
colnames(meanCount) <- c('Mock', 'flg22', 'flg22_SynCom33', 'flg22_SynCom35')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~cluster~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## scale
scaleCount <- meanCount %>%
  t %>%
  scale %>%
  t
scaleCount %<>% .[complete.cases(.), ]

## Cluster rows by Pearson correlation
hr <- scaleCount %>%
  t %>%
  cor(method = 'pearson') %>%
  {1 - .} %>%
  as.dist %>%
  hclust(method = 'complete')

## Clusters columns by Spearman correlation
hc <- scaleCount %>%
  cor(method = 'spearman') %>%
  {1 - .} %>%
  as.dist %>%
  hclust(method = 'complete')

cairo_pdf('heatmap.pdf')
heatmap.2(meanCount,
          Rowv = as.dendrogram(hr),
          Colv = as.dendrogram(hc),
          col = redgreen(100),
          scale = 'row',
          margins = c(7, 7),
          cexCol = 0.7,
          labRow = F,
          main = 'Heatmap.2',
          trace = 'none')
dev.off()

hc %>%
  as.dendrogram(method = 'average') %>%
  plot(main = 'Sample Clustering',
       ylab = 'Height')

hr %>%
  as.dendrogram(method = 'average') %>%
  plot(leaflab = 'none',
       main = 'Gene Clustering',
       ylab = 'Height')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~dynamic tree cut~~~~~~~~~~~~~~~~~~~~
clusDyn <- scaleCount %>%
  t %>%
  cor %>%
  {1 - .} %>%
  as.dist %T>%
  gc %>%
  as.matrix %T>%
  gc %>%
  cutreeDynamic(hr, distM = ., method = 'hybrid')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~cut trees by height ~~~~~~~~~~~~~~~~~~~~~~~~~
hclusth1.5 <- cutree(hr, h = 1.5)
hclusth1.0 <- cutree(hr, h = 1.0)
hclusth0.5 <- cutree(hr, h = 0.5)

cairo_pdf('genetree.pdf')
treeR <- hr %>%
  as.dendrogram(method = 'average')

plot(treeR,
     leaflab = 'none',
     main = 'Gene Clustering',
     ylab = 'Height')

cbind(hclusth0.5, hclusth1.0, hclusth1.5, clusDyn) %>%
  colored_bars(treeR,
               sort_by_labels_order = TRUE,
               y_shift = -0.1,
               rowLabels = c('h=0.5','h=1.0','h=1.5', 'Dynamic'),
               cex.rowLabels=0.7)

abline(h=1.5, lty = 2, col='grey')
abline(h=1.0, lty = 2, col='grey')
abline(h=0.5, lty = 2, col='grey')
dev.off()

cgenes <- c('AT1G14550.1', 'AT2G30750.1', 'AT2G19190.1')
hclusth1.5[cgenes]
hclusth1.0[cgenes]
hclusth0.5[cgenes]
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


##~~~~~~~~~~~~~~~~~~~~~~~~~~plot patterns~~~~~~~~~~~~~~~~~~~~~~~~
## join cluster and scaled normalized counts
clusterGene <- scaleCount %>%
  as.data.frame %>%
  rownames_to_column(var = 'ID') %>%
  as_tibble %>%
  {
    cl <- as.data.frame(hclusth1.5) %>%
      rownames_to_column(var = 'ID')
    inner_join(., cl)
  }

## plot core cluster
clusterCore <- clusterGene %>%
  group_by(hclusth1.5) %>%
  summarise_at(2:5, mean, na.rm = TRUE) %>% ## mean of each cluster
  mutate(hclusth1.5 = hclusth1.5 %>% paste0('cluster_', .)) %>%
  gather(Sample, NorExpress, Mock : flg22_SynCom35)
clusterCore$Sample %<>% factor(levels = c('Mock', 'flg22', 'flg22_SynCom33', 'flg22_SynCom35'), ordered = TRUE)

ggplot(clusterCore, aes(Sample, NorExpress, col = hclusth1.5, group = hclusth1.5)) +
  geom_point() +
  geom_line() +
  facet_wrap(. ~ hclusth1.5, ncol = 2) +
  ylab('Scaled counts') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave('hieracluster_1d5.jpg')
ggsave('hieracluster_1d5.pdf')

## plot all genes
clusterGenePlot <- clusterGene %>%
  gather(Sample, NorExpress, Mock : flg22_SynCom35) %>%
  mutate(hclusth1.5 = hclusth1.5 %>% paste0('cluster_', .))
clusterGenePlot$Sample %<>% factor(levels = c('Mock', 'flg22', 'flg22_SynCom33', 'flg22_SynCom35'), ordered = TRUE)

clusterCorePlot <- clusterCore %>% dplyr::mutate(ID = 1 : nrow(clusterCore))
ggplot(clusterGenePlot, aes(Sample, NorExpress, group = ID)) +
  geom_line(color = 'grey30', alpha = 0.01) +
  facet_wrap(. ~ hclusth1.5, ncol = 2) +
  geom_point(data = clusterCorePlot, aes(Sample, NorExpress, col = hclusth1.5, group = ID)) +
  geom_line(data = clusterCorePlot, aes(Sample, NorExpress, group = hclusth1.5, col = hclusth1.5)) +
  ylab('Scaled counts') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave('hieracluster_gene_1d5.pdf', width = 10, dpi = 320)
ggsave('hieracluster_gene_1d5.jpg', width = 10, dpi = 320)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~cluster cor phenotype~~~~~~~~~~~~~~~~~
traits <- data.frame(flg22 = c(0, 1, 1, 1),
                     SynCom33 = c(0, 0, 1, 0),
                     SynCom35 = c(0, 0, 0, 1),
                     rootlen = c(5.6, 1.1, 1.3, 4.9))

cores <- clusterGene %>%
  group_by(hclusth1.5) %>%
  summarise_at(2:5, mean, na.rm = TRUE) %>%
  mutate(hclusth1.5 = hclusth1.5 %>% paste0('cluster_', .)) %>%
  column_to_rownames(var = 'hclusth1.5') %>%
  t

moduleTraitCor <- cor(cores, traits, use = 'p')
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(traits))

traitPPlot <- moduleTraitPvalue %>%
  as.data.frame %>%
  rownames_to_column('cluster') %>%
  gather(trait, pvalue, flg22 : rootlen) %>%
  as_tibble

traitCorPlot <- moduleTraitCor %>%
  as.data.frame %>%
  rownames_to_column('cluster') %>%
  gather(trait, correlation, flg22 : rootlen) %>%
  as_tibble %>%
  mutate(x = rep(0 : (nrow(cores) - 1), each = ncol(cores))) %>%
  mutate(y = rep((ncol(cores) - 1) : 0, nrow(cores))) %>%
  inner_join(traitPPlot) %>%
  mutate(addtext = paste0(round(correlation, digit = 2),
                          '\n',
                          '(',
                          round(pvalue, digit = 2),
                          ')'))

ggplot(traitCorPlot, aes(x = x, y = y, fill = correlation)) +
  geom_tile() +
  scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name = 'RdYlBu')))(100),
                       breaks = seq(-1, 1, 0.5),
                       labels = format(seq(-1, 1, 0.5)),
                       limits = c(-1, 1)) +
  geom_text(aes(label = addtext)) +
  scale_x_continuous(breaks = 0 : 3, labels = c('flg22', 'SynCom33', 'SynCom35', 'rootlen')) +
  scale_y_continuous(breaks = 0 : 7, labels = paste0('cluster_', 8:1)) +
  xlab('Trait') +
  ylab('Cluster')
ggsave('hieracluster_1d5_trait.jpg')
ggsave('hieracluster_1d5_trait.pdf')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~heat map~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
scaleC <- rawCount %>%
  t %>%
  scale %>%
  t %>%
  as.data.frame %>%
  rownames_to_column('ID') %>%
  as_tibble %>%
  rename_at(-1, .funs = list(~paste0('Scale_', .)))

rawC <- rawCount %>%
  as.data.frame %>%
  rownames_to_column('ID') %>%
  as_tibble %>%
  rename_at(-1, .funs = list(~paste0('Raw_', .)))

degresC <- degres %>%
  select(ID, Flg22_vs_Mock_pvalue : Flg22_SynCom35_vs_Mock_log2FoldChange)

heatPlot <- rawC %>%
  inner_join(scaleC) %>%
  inner_join(degresC) %>%
  {
    cl <- as.data.frame(hclusth1.5) %>%
      rownames_to_column(var = 'ID')
    inner_join(., cl)
  } %>%
  slice(hr$order)

heatRawPlot <- heatPlot %>%
  select(ID, starts_with('Raw')) %>%
  gather(sample, raw, Raw_Mock_1 : Raw_Flg22_SynCom35_3) %>%
  mutate(x = rep(0 : 11, each = nrow(heatPlot))) %>%
  mutate(y = rep(0 : (nrow(heatPlot) - 1), 12))

heatScalePlot <- heatPlot %>%
  select(ID, starts_with('Scale')) %>%
  gather(sample, scale, Scale_Mock_1 : Scale_Flg22_SynCom35_3) %>%
  mutate(x = rep(0 : 11, each = nrow(heatPlot))) %>%
  mutate(y = rep(0 : (nrow(heatPlot) - 1), 12))

heatlog2FCPlot <- heatPlot %>%
  select(ID, ends_with('FoldChange')) %>%
  gather(sample, log2FC, Flg22_vs_Mock_log2FoldChange : Flg22_SynCom35_vs_Mock_log2FoldChange) %>%
  mutate(x = rep(0 : 2, each = nrow(heatPlot))) %>%
  mutate(y = rep(0 : (nrow(heatPlot) - 1), 3))

## sig |FC| > 1 and padj < 0.05
fcsig <- heatPlot %>%
  select(ends_with('FoldChange')) %>%
  abs %>%
  `>`(1)

padjsig <- heatPlot %>%
  select(ends_with('padj')) %>%
  abs %>%
  `<`(0.05) %>%
  apply(1:2, function(x) ifelse(is.na(x), FALSE, TRUE))

heatsigPlot <- (padjsig * fcsig) %>%
  as_tibble %>%
  gather(sample, sig, Flg22_vs_Mock_padj : Flg22_SynCom35_vs_Mock_padj) %>%
  mutate(x = rep(0 : 2, each = nrow(heatPlot))) %>%
  mutate(y = rep(0 : (nrow(heatPlot) - 1), 3))

heatGroupPlot <- heatPlot %>%
  select(ID, cluster = hclusth1.5) %>%
  mutate(x = 0) %>%
  mutate(y = 0 : (nrow(heatPlot) - 1))

theme_flg22 <- function(...) {
  theme_bw() %+replace%
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks.length = unit(0, 'mm'),
          axis.line = element_blank(),
          panel.spacing = unit(0, 'mm'),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), 'line'),
          legend.spacing = unit(0, 'mm'),
          ...)
}

ggplot(heatRawPlot, aes(x = x, y = y, fill = log2(raw))) +
  geom_tile() +
  scale_fill_gradientn(colours = colorRampPalette(brewer.pal(n = 7, name = 'GnBu'))(100), name = 'log2(count)') +
  scale_x_continuous(breaks = 0 : 11,
                     labels = rep(c('Mock', 'flg22', 'flg22_SynCom33', 'flg22_SynCom35'), each = 3) %>%
                       paste(rep(1 : 3, 4), sep = '_')) +
  theme_flg22(legend.position = 'top',
              axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(heatScalePlot, aes(x = x, y = y, fill = scale)) +
  geom_tile() +
  scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name = 'RdYlBu')))(100), name = 'scale(count)') +
  scale_x_continuous(breaks = 0 : 11,
                     labels = rep(c('Mock', 'flg22', 'flg22_SynCom33', 'flg22_SynCom35'), each = 3) %>%
                       paste(rep(1 : 3, 4), sep = '_')) +
  theme_flg22(legend.position = 'top',
              axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(heatlog2FCPlot, aes(x = x, y = y, fill = log2FC)) +
  geom_tile() +
  scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 8, name = 'PiYG')))(100), name = 'log2(FoldChange)') +
  scale_x_continuous(breaks = 0 : 2,
                     labels = paste(c('flg22', 'flg22_SynCom33', 'flg22_SynCom35'), 'vs. Mock')) +
  theme_flg22(legend.position = 'top',
              axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(heatsigPlot, aes(x = x, y = y)) +
  geom_tile(aes(fill = factor(sig))) +
  scale_fill_manual(name = 'significant', labels = c('No', 'Yes'), values = c('grey', 'red'))

ggplot(heatGroupPlot, aes(x = x, y = y)) +
  geom_tile(aes(fill = factor(cluster))) +
  scale_x_continuous(breaks = 0,
                     labels = 'group') +
  theme_flg22(legend.position = 'left',
              axis.text.x = element_text(angle = 90, hjust = 1))



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################################################################
