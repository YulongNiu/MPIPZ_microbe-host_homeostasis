###########################replot heatmap##############################
library('tidyverse')
library('ComplexHeatmap')
library('limma')
library('DESeq2')
library('RColorBrewer')

setwd('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/')

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~load WER and Col0 data~~~~~~~~~~~~~~~~~~~
load('degres_condi_Mock.RData')
kmeansRes <- read_csv('kmeans10.csv') %>%
  select(ID, cl)

scaleCWER <- rldData %>%
  as.data.frame %>%
  rownames_to_column('ID') %>%
  as_tibble %>%
  select(contains('_')) %>%
  t %>%
  scale %>%
  t %>%
  as_tibble %>%
  bind_cols(ID = rownames(rldData)) %>%
  inner_join(kmeansRes)

load('degres_condi_Mock_1stadd.RData')
kmeansRes <- read_csv('kmeans10_1stadd.csv') %>%
  select(ID, cl)

scaleCCol0 <- rldData %>%
  as.data.frame %>%
  rownames_to_column('ID') %>%
  as_tibble %>%
  select(contains('_')) %>%
  t %>%
  scale %>%
  t %>%
  as_tibble %>%
  bind_cols(ID = rownames(rldData)) %>%
  inner_join(kmeansRes) %>%
  dplyr::select(matches('Mock_\\d|HKSynCom33_\\d|HKSynCom35_\\d'), matches('Mock_Flg22_\\d|HKSynCom33_Flg22_\\d|HKSynCom35_Flg22_\\d'), everything())
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~combine heatmap~~~~~~~~~~~~~~~~~~~~~~~~~~~
sigCol0 <- read_csv('kmeans10_1stadd_sig.csv') %>%
  select(ID) %>%
  inner_join(scaleCCol0 %>% select(ID)) %>%
  inner_join(scaleCWER %>% select(ID))

scaleCCol0sig <- scaleCCol0 %>%
  inner_join(sigCol0)

scaleCWERSig <- scaleCWER %>%
  inner_join(sigCol0)

ht_list <- Heatmap(matrix = scaleCCol0sig %>% select(contains('_')),
        name = 'Scaled Counts',
        row_order = order(scaleCCol0sig$cl) %>% rev,
        row_split = scaleCCol0sig$cl,
        row_gap = unit(2, "mm"),
        column_order = 1 : 40,
        column_split = rep(c('Mock/HKSynCom', 'Non-sup', 'Sup'), c(24, 8, 8)),
        show_column_names = FALSE,
        col = colorRampPalette(rev(brewer.pal(n = 10, name = 'Spectral'))[c(-3, -4, -6, -7)])(100),
        column_title_gp = gpar(fontsize = 10),
        use_raster = FALSE) +
  Heatmap(matrix = scaleCWERSig %>% select(contains('_')),
          name = 'Scaled Counts',
          column_order = 1 : 12,
          column_split = rep(c('Mock', 'Mock+flg22', 'Non-sup+flg22', 'Sup+flg22'), each = 3),
          show_column_names = FALSE,
          col = colorRampPalette(rev(brewer.pal(n = 10, name = 'Spectral'))[c(-3, -4, -6, -7)])(100),
          column_title_gp = gpar(fontsize = 6),
          use_raster = FALSE)

pdf('kmeans10_merge.pdf', width = 10)
draw(ht_list)
dev.off()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#######################################################################
