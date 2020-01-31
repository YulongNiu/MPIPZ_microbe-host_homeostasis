###########################replot heatmap##############################
library('tidyverse')
library('ComplexHeatmap')
library('limma')
library('DESeq2')
library('RColorBrewer')

setwd('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/')

load('degres_condi_Mock_1stadd.RData')

##~~~~~~~~~~~~~~~~~~~~~~useful funcs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
meanFlg22 <- function(v) {

  require('magrittr')

  res <- v %>%
    split(rep(1 : 10, each = 4)) %>%
    sapply(mean, na.rm = TRUE)

  return(res)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~select DEGs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
wholeDEG <- read_csv('eachGroup_vs_Mock_k_1stadd.csv')
kmeansRes <- read_csv('kmeans10_1stadd.csv') %>%
  select(ID, cl)

fcsig <- wholeDEG %>%
  select(ends_with('FoldChange')) %>%
  transmute_all(list(~ case_when(. > log2(1.5) ~ 1,
                                 . < -log2(1.5) ~ -1,
                                 TRUE ~ 0)))
padjsig <- wholeDEG %>%
  select(ends_with('padj')) %>%
  `<`(0.05) %>%
  as_tibble %>%
  transmute_all(list(~ if_else(is.na(.), FALSE, .)))

heatsig <- (padjsig * fcsig) %>%
  as_tibble %>%
  abs %>%
  rowSums %>%
  {. >= 1} %>%
  which %>%
  dplyr::slice(wholeDEG, .) %>%
  inner_join(kmeansRes)

heatsig %>%
  mutate_at(c('Gene', 'Description'), .funs = list(~if_else(is.na(.), '', .))) %>%
  write_csv('kmeans10_1stadd_sig.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~heatmap~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## rlog transformed
rawC <- rldData %>%
  as.data.frame %>%
  rownames_to_column('ID') %>%
  as_tibble %>%
  dplyr::select(matches('Mock_\\d|HKSynCom33_\\d|HKSynCom35_\\d'), matches('Mock_Flg22_\\d|HKSynCom33_Flg22_\\d|HKSynCom35_Flg22_\\d'), everything()) %>%
  inner_join(heatsig %>% select(ID, cl))
  ## inner_join(kmeansRes) ## all transcripts

## scale counts
scaleC <- rawC %>%
  select(contains('_')) %>%
  t %>%
  scale %>%
  t %>%
  as_tibble %>%
  bind_cols(rawC %>% select(ID, cl))

## 1. sample annotation
col_flg22 <- HeatmapAnnotation(Flg22 = c(rep(c('No', 'Yes'), each = 12),
                                     rep(c('No', 'Yes', 'No', 'Yes'), each = 4)),
                          col = list(Flg22 = c('Yes' = 'grey', 'No' = 'white')),
                          gp = gpar(col = 'black'))

col_syncom <- HeatmapAnnotation(SynCom = c(rep(c('Mock', 'HK', 'Mock', 'HK'), c(4, 8, 4, 8)),
                                       rep('Live', 16)),
                            col = list(SynCom = c('Mock' = 'grey80', 'HK' = 'white', 'Live' = 'grey50')),
                            gp = gpar(col = 'black'))

## 2. DEG annotation
sigMat <- (padjsig * fcsig) %>%
  as_tibble %>%
  setNames(names(.) %>% substr(., start = 1, stop = nchar(.) - 5)) %>%
  mutate(ID = wholeDEG$ID) %>%
  inner_join(scaleC %>% select(ID), .) %T>%
  {(sum(.$ID == scaleC$ID) == nrow(.)) %>% print} %>%
  transmute_at(.var = vars(contains('vs')),
               list(~ case_when(. == -1 ~ 'down',
                                 . == 0 ~'no',
                                . == 1 ~ 'up'))) %>%
  as.matrix

## flg22 matrix
flg22Mat <- sigMat[, 1:5]
colnames(flg22Mat) <- c('Mock+flg22 vs. Mock',
                        'HKSynCom33+flg22 vs. HKSynCom33',
                        'HKSynCom35+flg22 vs. HKSynCom35',
                        'SynCom33+flg22 vs. SynCom33',
                        'SynCom35+flg22 vs. SynCom35')

bacMat <- sigMat[, 6:11]
colnames(bacMat) <- c('HKSynCom33 vs. Mock',
                      'HKSynCom35 vs. Mock',
                      'SynCom33 vs. Mock',
                      'SynCom35 vs. Mock',
                      'SynCom33 vs. SynCom35',
                      'HKSynCom33 vs. HkSynCom35')

mixMat <- sigMat[, 12:21]
colnames(mixMat) <- c('HKSynCom33+flg22 vs. Mock',
                      'HKSynCom33+flg22 vs. Mock',
                      'SynCom33+flg22 vs. Mock',
                      'SynCom35+flg22 vs. Mock',
                      'SynCom33+flg22 vs. SynCom35+flg22',
                      'HKSynCom33+flg22 vs. HKSynCom35+flg22',
                      'HKSynCom35+flg22 vs. Mock+flg22',
                      'HKSynCom33+flg22 vs. Mock+flg22',
                      'SynCom33+flg22 vs. Mock+flg22',
                      'SynCom35+flg22 vs. Mock+flg22')

hkMat <- sigMat[, 22:23]
colnames(hkMat) <- c('SynCom33 vs. HKSynCom33',
                     'SynCom35 vs. HKSynCom35')

ht_list <- Heatmap(matrix = scaleC %>% select(contains('_')),
                   name = 'Scaled Counts',
                   row_order = order(scaleC$cl) %>% rev,
                   row_split = scaleC$cl,
                   row_gap = unit(2, "mm"),
                   column_order = 1 : 40,
                   column_split = rep(c('Mock/HKSynCom', 'Non-sup', 'Sup'), c(24, 8, 8)),
                   show_column_names = FALSE,
                   col = colorRampPalette(rev(brewer.pal(n = 10, name = 'Spectral'))[c(-3, -4, -6, -7)])(100),
                   top_annotation = c(col_flg22, col_syncom),
                   use_raster = FALSE) +
  Heatmap(flg22Mat,
          col = c('down' = 'blue', 'no' = 'white', 'up' = 'red'),
          column_names_gp = gpar(fontsize = 5),
          heatmap_legend_param = list(title = 'DEGs'),
          cluster_columns = FALSE,
          use_raster = FALSE) +
  Heatmap(bacMat,
        col = c('down' = 'blue', 'no' = 'white', 'up' = 'red'),
        column_names_gp = gpar(fontsize = 5),
        cluster_columns = FALSE,
        use_raster = FALSE,
        show_heatmap_legend = FALSE) +
  Heatmap(mixMat,
          col = c('down' = 'blue', 'no' = 'white', 'up' = 'red'),
          column_names_gp = gpar(fontsize = 5),
          cluster_columns = FALSE,
          use_raster = FALSE,
          show_heatmap_legend = FALSE) +
  Heatmap(hkMat,
          col = c('down' = 'blue', 'no' = 'white', 'up' = 'red'),
          column_names_gp = gpar(fontsize = 5),
          cluster_columns = FALSE,
          use_raster = FALSE,
          show_heatmap_legend = FALSE)


pdf('kmeans10_heatmap_1stadd_sig2.pdf')
draw(ht_list)
dev.off()
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#######################################################################
