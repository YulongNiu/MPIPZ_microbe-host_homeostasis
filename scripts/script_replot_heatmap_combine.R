##########################WER and Col0##############################
library('tidyverse')
library('ComplexHeatmap')
library('limma')
library('DESeq2')
library('RColorBrewer')

setwd('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/removeZero')

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
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~Col0 arga DEGs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

sigMat <- (padjsig * fcsig) %>%
  as_tibble %>%
  setNames(names(.) %>% substr(., start = 1, stop = nchar(.) - 5)) %>%
  mutate(ID = wholeDEG$ID) %>%
  inner_join(scaleCCol0 %>% select(ID), .) %T>%
  {(sum(.$ID == scaleCCol0$ID) == nrow(.)) %>% print} %>%
  transmute_at(.var = vars(contains('vs')),
               list(~ case_when(. == -1 ~ 'down',
                                . == 0 ~'no',
                                . == 1 ~ 'up'))) %>%
  mutate(ID = scaleCCol0$ID)

## flg22 matrix
flg22Col <- sigMat %>% select(1:3, ID)
colnames(flg22Col)[1:3] <- c('Mock+flg22 vs. Mock',
                             'HKSynCom33+flg22 vs. HKSynCom33',
                             'HKSynCom35+flg22 vs. HKSynCom35')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~Iron bac response~~~~~~~~~~~~~~~~~~~~~~~~
ironBacSig <- read_csv('/extDisk1/RESEARCH/MPIPZ_CJ_RNASeq/results/eachGroup_mergeDay8_deg_sig.csv') %>%
  select(ID, cl) %>%
  filter(cl %in% c(10, 3))

ironRespSig <- read_csv('/extDisk1/RESEARCH/MPIPZ_CJ_RNASeq/results/eachGroup_mergeDay8_deg_sig.csv') %>%
  select(ID, cl) %>%
  filter(cl %in% c(4, 8))
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~combine heatmap~~~~~~~~~~~~~~~~~~~~~~~~~~~
sigCol0 <- read_csv('kmeans10_1stadd_sig.csv') %>%
  select(ID) %>%
  inner_join(scaleCCol0 %>% select(ID)) %>%
  inner_join(scaleCWER %>% select(ID))

scaleCCol0sig <- scaleCCol0 %>%
  inner_join(sigCol0)

scaleCWERSig <- scaleCWER %>%
  inner_join(sigCol0)

flg22Col0Sig <- flg22Col %>%
  inner_join(sigCol0)

ironCol0Sig <- sigCol0 %>%
  left_join(., ironBacSig) %>%
  dplyr::rename(ironbac = cl) %>%
  mutate(ironbac = case_when(ironbac == 10 ~ 'bacup',
                             ironbac == 3 ~ 'bacdown',
                             is.na(ironbac) ~ 'bacno')) %>%
  left_join(., ironRespSig) %>%
  dplyr::rename(ironresp = cl) %>%
  mutate(ironresp = case_when(ironresp == 8 ~ 'respup',
                              ironresp == 4 ~ 'respdown',
                              is.na(ironresp) ~ 'respno')) %>%
  inner_join(scaleCCol0sig %>% select(ID), .)

all.equal(sum(scaleCCol0sig$ID == scaleCWERSig$ID),
          sum(scaleCCol0sig$ID == flg22Col0Sig$ID),
          sum(scaleCCol0sig$ID == ironBacCol0Sig$ID),
          nrow(scaleCCol0sig))

ht_list <- Heatmap(matrix = scaleCCol0sig %>% select(contains('_')),
        name = 'Scaled Counts',
        ## row_order = order(scaleCCol0sig$cl) %>% rev,
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
          use_raster = FALSE) +
  Heatmap(flg22Col0Sig %>% select(-ID),
          col = c('down' = 'blue', 'no' = 'white', 'up' = 'red'),
          column_names_gp = gpar(fontsize = 5),
          heatmap_legend_param = list(title = 'DEGs'),
          cluster_columns = FALSE,
          use_raster = FALSE) +
  Heatmap(ironCol0Sig %>% select(ironbac),
          col = c('bacup' = 'purple', 'bacno' = 'white', 'bacdown' = 'green3'),
          column_names_gp = gpar(fontsize = 5),
          heatmap_legend_param = list(title = 'IronBac'),
          cluster_columns = FALSE,
          use_raster = FALSE) +
  Heatmap(ironCol0Sig %>% select(ironresp),
          col = c('respup' = 'purple', 'respno' = 'white', 'respdown' = 'green3'),
          column_names_gp = gpar(fontsize = 5),
          heatmap_legend_param = list(title = 'IronResp'),
          cluster_columns = FALSE,
          use_raster = FALSE)

## filePrefix <- 'kmeans10_heatmap_WER_Col02'
filePrefix <- 'kmeans10_heatmap_WER_Col02_Iron2'

pdf(paste0(filePrefix, '.pdf'))
draw(ht_list)
dev.off()

system(paste0('convert -density 1200 ', paste0(filePrefix, '.pdf'), ' ', paste0(filePrefix, '.jpg')))
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#######################################################################

##########################Col0 agar soil##############################
library('tidyverse')
library('ComplexHeatmap')
library('limma')
library('DESeq2')
library('RColorBrewer')

setwd('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/removeZero')

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~load WER and Col0 data~~~~~~~~~~~~~~~~~~~
load('degres_condi_Mock_soil.RData')

scaleCSoil <- rldData %>%
  as.data.frame %>%
  rownames_to_column('ID') %>%
  as_tibble %>%
  select(contains('_')) %>%
  t %>%
  scale %>%
  t %>%
  as_tibble %>%
  bind_cols(ID = rownames(rldData))

load('degres_condi_Mock_1stadd.RData')
kmeansRes <- read_csv('kmeans10_1stadd.csv') %>%
  select(ID, cl)

scaleCAgar <- rldData %>%
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
  inner_join(scaleCAgar %>% select(ID)) %>%
  inner_join(scaleCSoil %>% select(ID))

scaleCAgarsig <- scaleCAgar %>%
  inner_join(sigCol0)

scaleCSoilSig <- scaleCSoil %>%
  inner_join(sigCol0)

ht_list <- Heatmap(matrix = scaleCAgarsig %>% select(contains('_')),
        name = 'Scaled Counts',
        row_order = order(scaleCAgarsig$cl) %>% rev,
        row_split = scaleCAgarsig$cl,
        row_gap = unit(2, "mm"),
        column_order = 1 : 40,
        column_split = rep(c('Mock/HKSynCom', 'Non-sup', 'Sup'), c(24, 8, 8)),
        show_column_names = FALSE,
        col = colorRampPalette(rev(brewer.pal(n = 10, name = 'Spectral'))[c(-3, -4, -6, -7)])(100),
        column_title_gp = gpar(fontsize = 10),
        use_raster = FALSE) +
  Heatmap(matrix = scaleCSoilSig %>% select(contains('_')),
          name = 'Scaled Counts',
          column_order = 1 : 9,
          column_split = rep(c('Mock', 'Non-sup', 'Sup'), each = 3),
          show_column_names = FALSE,
          col = colorRampPalette(rev(brewer.pal(n = 10, name = 'Spectral'))[c(-3, -4, -6, -7)])(100),
          column_title_gp = gpar(fontsize = 6),
          use_raster = FALSE)

filePrefix <- 'kmeans10_heatmap_soil_agar'

pdf(paste0(filePrefix, '.pdf'))
draw(ht_list)
dev.off()

system(paste0('convert -density 1200 ', paste0(filePrefix, '.pdf'), ' ', paste0(filePrefix, '.jpg')))
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#######################################################################
