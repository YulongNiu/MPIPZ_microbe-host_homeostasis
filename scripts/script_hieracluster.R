######################hierarchical clustering####################
setwd('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/')

load('degres_condi_Mock.RData')

library('readr')
library('dplyr')
library('magrittr')
library('tibble')
library('gplots')
library('dendextend')
library('dynamicTreeCut')
library('ggplot2')
library('tidyr')
library('DESeq2')

##~~~~~~~~~~~~~~~~~~~~~~useful funcs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
meanFlg22 <- function(v) {

  require('magrittr')

  res <- v %>%
    split(rep(1 : 4, each = 3)) %>%
    sapply(mean, na.rm = TRUE)

  return(res)
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

cbind(hclusth0.5, hclusth1.0, hclusth1.5) %>%
  colored_bars(treeR,
               sort_by_labels_order = TRUE,
               y_shift = -0.1,
               rowLabels = c('h=0.5','h=1.0','h=1.5'),
               cex.rowLabels=0.7)

abline(h=1.5, lty = 2, col='grey')
abline(h=1.0, lty = 2, col='grey')
abline(h=0.5, lty = 2, col='grey')
dev.off()

anno %>%
  mutate(cluster = hclusth0.5) %>%
  write_csv('tmp1.csv')
cgenes <- c('AT1G14550.1', 'AT2G30750.1', 'AT2G19190.1')
hclusth1.5[cgenes]
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~dynamic tree cut~~~~~~~~~~~~~~~~~~~~
clusDyn <- scaleCount %>%
  t %>%
  cor %>%
  {1 - .} %>%
  as.dist

cutreeDynamic(hr, clusDyn, method = 'hybrid')
clusDyn <- cutreeDynamic(hr, distM = as.matrix(as.dist(1-cor(t(scaledata)))), method = "hybrid")
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


##~~~~~~~~~~~~~~~~~~~~~~~~~~plot patterns~~~~~~~~~~~~~~~~~~~~~~~~
## plot cores
clusterCore <- scaleCount %>%
  as.data.frame %>%
  rownames_to_column(var = 'ID') %>%
  as_tibble %>%
  {
    cl <- as.data.frame(hclusth1.5) %>%
      rownames_to_column(var = 'ID')
    inner_join(., cl)
  } %>% ## join cluster and scaled normalized counts
  group_by(hclusth1.5) %>%
  summarise_at(2:5, mean, na.rm = TRUE) %>% ## mean of each cluster
  mutate(hclusth1.5 = hclusth1.5 %>% paste0('cluster_', .)) %>%
  gather(Sample, NorExpress, Mock : flg22_SynCom35)

clusterCore$Sample %<>% factor(levels = c('Mock', 'flg22', 'flg22_SynCom33', 'flg22_SynCom35'), ordered = TRUE)

ggplot(clusterCore, aes(Sample, NorExpress, col = hclusth1.5, group = hclusth1.5)) +
  geom_point() +
  geom_line() +
  ylab('Scaled counts') +
  facet_wrap(. ~ hclusth1.5, ncol = 2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave('hieracluster.pdf')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################################################################
