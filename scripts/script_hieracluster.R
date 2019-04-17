######################hierarchical clustering####################

setwd('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/')

library('readr')
library('dplyr')
library('magrittr')
library('tibble')
library('gplots')

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
rawCount <- read_csv('eachGroup_vs_Mock_k.csv',
                     col_types = cols(Chromosome = col_character())) %>%
  select(ID, Mock_1 : Flg22_SynCom35_3) %>%
  column_to_rownames('ID') %>%
  as.matrix

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
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################################################################
