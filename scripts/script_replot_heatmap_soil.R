
## originally by Yulong Niu
## yulong.niu@hotmail.com

###########################replot heatmap##############################
library('tidyverse')
library('ComplexHeatmap')
library('limma')
library('DESeq2')
library('RColorBrewer')

setwd('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/removeZero')

load('degres_condi_Mock_soil.RData')

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~select DEGs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
wholeDEG <- read_csv('eachGroup_vs_Mock_k_soil.csv')
kmeansRes <- read_csv('kmeans10_soil.csv') %>%
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
  write_csv('kmeans10_soil_sig.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~heatmap~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## rlog transformed
rawC <- rldData %>%
  as.data.frame %>%
  rownames_to_column('ID') %>%
  as_tibble %>%
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

## 1. DEG annotation
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

colnames(sigMat) <- c('SynCom33 vs. Mock',
                      'SynCom35 vs. Mock',
                      'SynCom33 vs. SynCom35')

ht_list <- Heatmap(matrix = scaleC %>% select(contains('_')),
                   name = 'Scaled Counts',
                   row_order = order(scaleC$cl) %>% rev,
                   row_split = scaleC$cl,
                   row_gap = unit(2, "mm"),
                   column_order = 1 : 9,
                   column_split = rep(c('Mock', 'Non-sup', 'Sup'), each = 3),
                   show_column_names = FALSE,
                   col = colorRampPalette(rev(brewer.pal(n = 10, name = 'Spectral'))[c(-3, -4, -6, -7)])(100),
                   column_title_gp = gpar(fontsize = 7),
                   use_raster = FALSE) +
  Heatmap(sigMat,
          col = c('down' = 'blue', 'no' = 'white', 'up' = 'red'),
          column_names_gp = gpar(fontsize = 5),
          heatmap_legend_param = list(title = 'DEGs'),
          cluster_columns = FALSE,
          use_raster = FALSE)

filePrefix <- 'kmeans10_heatmap_soil_sig_DEG'

pdf(paste0(filePrefix, '.pdf'))
draw(ht_list)
dev.off()

system(paste0('convert -density 1200 ', paste0(filePrefix, '.pdf'), ' ', paste0(filePrefix, '.jpg')))
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#######################################################################


