
## originally by Yulong Niu
## yulong.niu@hotmail.com

##########################WER and Col0##############################
library('clusterProfiler')
library('tidyverse')
library('ComplexHeatmap')
library('limma')
library('DESeq2')
library('RColorBrewer')
library('foreach')
library('magrittr')
library('org.At.tair.db')

basepath <- '/extDisk1/RESEARCH/'

basepath %>%
  file.path('MPIPZ_KaWai_RNASeq/results/removeZero') %>%
  setwd

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~load WER and Col0 data~~~~~~~~~~~~~~~~~~~
load('degres_condi_Mock.RData')
kmeansRes <- read_csv('kmeans10.csv') %>%
  dplyr::select(ID, cl)

scaleCWER <- rldData %>%
  as.data.frame %>%
  rownames_to_column('ID') %>%
  as_tibble %>%
  dplyr::select(contains('_')) %>%
  t %>%
  scale %>%
  t %>%
  as_tibble %>%
  bind_cols(ID = rownames(rldData)) %>%
  inner_join(kmeansRes)

load('degres_condi_Mock_1stadd.RData')
kmeansRes <- read_csv('kmeans10_1stadd.csv') %>%
  dplyr::select(ID, cl)

scaleCCol0 <- rldData %>%
  as.data.frame %>%
  rownames_to_column('ID') %>%
  as_tibble %>%
  dplyr::select(contains('_')) %>%
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
  dplyr::select(ID, cl)

fcsig <- wholeDEG %>%
  dplyr::select(ends_with('FoldChange')) %>%
  transmute_all(list(~ case_when(. > log2(1.5) ~ 1,
                                 . < -log2(1.5) ~ -1,
                                 TRUE ~ 0)))
padjsig <- wholeDEG %>%
  dplyr::select(ends_with('padj')) %>%
  `<`(0.05) %>%
  as_tibble %>%
  transmute_all(list(~ if_else(is.na(.), FALSE, .)))

sigMat <- (padjsig * fcsig) %>%
  as_tibble %>%
  setNames(names(.) %>% substr(., start = 1, stop = nchar(.) - 5)) %>%
  mutate(ID = wholeDEG$ID) %>%
  inner_join(scaleCCol0 %>% dplyr::select(ID), .) %T>%
  {(sum(.$ID == scaleCCol0$ID) == nrow(.)) %>% print} %>%
  transmute_at(.var = vars(contains('vs')),
               list(~ case_when(. == -1 ~ 'down',
                                . == 0 ~'no',
                                . == 1 ~ 'up'))) %>%
  mutate(ID = scaleCCol0$ID)

## flg22 matrix
flg22Col <- sigMat %>% dplyr::select(1:5, ID)
colnames(flg22Col)[1:5] <- c('Mock+flg22 vs. Mock',
                             'HKSynCom33+flg22 vs. HKSynCom33',
                             'HKSynCom35+flg22 vs. HKSynCom35',
                             'SynCom33+flg22 vs. SynCom33',
                             'SynCom35+flg22 vs. SynCom35')

hkCol <- sigMat %>% dplyr::select(22:23, ID)
colnames(hkCol)[1:2] <- c('SynCom33 vs. HKSynCom33',
                          'SynCom35 vs. HKSynCom35')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~Iron bac response~~~~~~~~~~~~~~~~~~~~~~~~
## ironBacSig <- read_csv('/extDisk1/RESEARCH/MPIPZ_CJ_RNASeq/results/eachGroup_mergeDay8_deg_sig.csv') %>%
##   select(ID, cl) %>%
##   filter(cl %in% c(10, 3))

ironHKLive <- basepath %>%
  file.path('MPIPZ_CJ_RNASeq/results/eachGroup_mergeDay8.csv') %>%
  read_csv

fcsigIron <- ironHKLive %>%
  dplyr::select(ends_with('FoldChange')) %>%
  transmute_all(list(~ case_when(. > log2(1.5) ~ 1,
                                 . < -log2(1.5) ~ -1,
                                 TRUE ~ 0)))
padjsigIron <- ironHKLive %>%
  dplyr::select(ends_with('padj')) %>%
  `<`(0.05) %>%
  as_tibble %>%
  transmute_all(list(~ if_else(is.na(.), FALSE, .)))

sigMatIronDay8 <- (padjsigIron * fcsigIron) %>%
  as_tibble %>%
  setNames(names(.) %>% substr(., start = 1, stop = nchar(.) - 5)) %>%
  transmute_at(.var = vars(contains('vs')),
               list(~ case_when(. == -1 ~ 'bacdown',
                                . == 0 ~'bacno',
                                . == 1 ~ 'bacup'))) %>%
  mutate(ID = ironHKLive$ID)
colnames(sigMatIronDay8)[1:8] %<>% paste0('_Day8')


ironHKLive <- basepath %>%
  file.path('MPIPZ_CJ_RNASeq/results/eachGroup_mergeDay15.csv') %>%
  read_csv

fcsigIron <- ironHKLive %>%
  dplyr::select(ends_with('FoldChange')) %>%
  transmute_all(list(~ case_when(. > log2(1.5) ~ 1,
                                 . < -log2(1.5) ~ -1,
                                 TRUE ~ 0)))
padjsigIron <- ironHKLive %>%
  dplyr::select(ends_with('padj')) %>%
  `<`(0.05) %>%
  as_tibble %>%
  transmute_all(list(~ if_else(is.na(.), FALSE, .)))

sigMatIronDay15 <- (padjsigIron * fcsigIron) %>%
  as_tibble %>%
  setNames(names(.) %>% substr(., start = 1, stop = nchar(.) - 5)) %>%
  transmute_at(.var = vars(contains('vs')),
               list(~ case_when(. == -1 ~ 'bacdown',
                                . == 0 ~'bacno',
                                . == 1 ~ 'bacup'))) %>%
  mutate(ID = ironHKLive$ID)
colnames(sigMatIronDay15)[1:8] %<>% paste0('_Day15')


ironRespSig <- basepath %>%
  file.path('MPIPZ_CJ_RNASeq/results/eachGroup_mergeDay8_deg_sig.csv') %>%
  read_csv %>%
  dplyr::select(ID, cl) %>%
  filter(cl %in% c(4, 8))
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~Castrillo flg22 longtime~~~~~~~~~~~~~~~~
CastrilloSig <- basepath %>%
  file.path('MPIPZ_KaWai_RNASeq/flg22_crossref/Castrillo_2017/DEGsflg22_singleend.csv') %>%
  read_csv %>%
  dplyr::select(ID, flg22_vs_Col0_log2FoldChange) %>%
  dplyr::rename(Castrillo_log2FC = flg22_vs_Col0_log2FoldChange)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~Volz flg22 shorttime~~~~~~~~~~~~~~~~
VolzSig <- basepath %>%
  file.path('MPIPZ_KaWai_RNASeq/flg22_crossref/Volz_2019/DEGs_pairend.csv') %>%
  read_csv %>%
  dplyr::select(ID, Col0_vs_flg22_log2FoldChange) %>%
  dplyr::rename(Volz_log2FC = Col0_vs_flg22_log2FoldChange) %>%
  mutate(Volz_log2FC = -Volz_log2FC)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~Paulo flg22~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PauloSig <- read_delim('Paulo_RNASeq_flg22.txt', delim = '\t')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~Paulo HK vs living~~~~~~~~~~~~~~~~~~~~~~~
PauloHKSig <- read_csv('Paulo_RNASeq_hklive_noflg22.csv') %>%
  dplyr::filter(abs(Fold) > 1.5, FDR < 0.05) %>%
  dplyr::select(Gene, logFC)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~Stringlis flg22~~~~~~~~~~~~~~~~~~~~~~~~~~~
StringlisDEGs <- basepath %>%
  file.path('MPIPZ_KaWai_RNASeq/flg22_crossref/Stringlis_2018/results/eachGroup_DEGs_Stringlis.csv') %>%
  read_csv

fcsigStringlis <- StringlisDEGs %>%
  dplyr::select(ends_with('FoldChange')) %>%
  transmute_all(list(~ case_when(. > log2(1.5) ~ 1,
                                 . < -log2(1.5) ~ -1,
                                 TRUE ~ 0)))

padjsigStringlis <- StringlisDEGs %>%
  dplyr::select(ends_with('padj')) %>%
  `<`(0.05) %>%
  as_tibble %>%
  transmute_all(list(~ if_else(is.na(.), FALSE, .)))

sigMatStringlis <- (padjsigStringlis * fcsigStringlis) %>%
  as_tibble %>%
  setNames(names(.) %>% substr(., start = 1, stop = nchar(.) - 5)) %>% {
    flg22 <- transmute_at(., .var = vars(starts_with('flg22')),
                       list(~ case_when(. == -1 ~ 'down',
                                        . == 0 ~'no',
                                        . == 1 ~ 'up')))

    bac <- transmute_at(., .var = vars(starts_with('WCS417')),
                     list(~ case_when(. == -1 ~ 'bacdown',
                                      . == 0 ~'bacno',
                                      . == 1 ~ 'bacup')))

    return(bind_cols(flg22, bac))
  } %>%
mutate(ID = StringlisDEGs$ID)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~combine heatmap~~~~~~~~~~~~~~~~~~~~~~~~~~~
sigCol0 <- read_csv('kmeans10_1stadd_sig.csv') %>%
  dplyr::select(ID) %>%
  inner_join(scaleCCol0 %>% dplyr::select(ID)) %>%
  inner_join(scaleCWER %>% dplyr::select(ID)) %>%
  mutate(Gene = ID %>% substring(first = 1, last = nchar(.) - 2))

scaleCCol0sig <- scaleCCol0 %>%
  inner_join(sigCol0)

scaleCWERSig <- scaleCWER %>%
  inner_join(sigCol0)

flg22Col0Sig <- flg22Col %>%
  inner_join(sigCol0)

hkCol0Sig <- hkCol %>%
  inner_join(sigCol0)

ironHKLiveSigDay8 <- sigCol0 %>%
  left_join(., sigMatIronDay8) %>%
  mutate_all(.funs = list(~if_else(is.na(.), 'bacno', .))) %>%
  inner_join(scaleCCol0sig %>% dplyr::select(ID), .)

ironHKLiveSigDay15 <- sigCol0 %>%
  left_join(., sigMatIronDay15) %>%
  mutate_all(.funs = list(~if_else(is.na(.), 'bacno', .))) %>%
  inner_join(scaleCCol0sig %>% dplyr::select(ID), .)

ironCol0Sig <- sigCol0 %>%
  ## left_join(., ironBacSig) %>%
  ## dplyr::rename(ironbac = cl) %>%
  ## mutate(ironbac = case_when(ironbac == 10 ~ 'bacup',
  ##                            ironbac == 3 ~ 'bacdown',
  ##                            is.na(ironbac) ~ 'bacno')) %>%
  left_join(., ironRespSig) %>%
  dplyr::rename(ironresp = cl) %>%
  mutate(ironresp = case_when(ironresp == 8 ~ 'respup',
                              ironresp == 4 ~ 'respdown',
                              is.na(ironresp) ~ 'respno')) %>%
  inner_join(scaleCCol0sig %>% dplyr::select(ID), .)

CastrilloCol0Sig <- sigCol0 %>%
  left_join(., CastrilloSig) %>%
  mutate(Castrillo_log2FC = case_when(Castrillo_log2FC > 0 ~ 'up',
                                      Castrillo_log2FC < 0 ~ 'down',
                                      is.na(Castrillo_log2FC) ~ 'no')) %>%
  inner_join(scaleCCol0sig %>% dplyr::select(ID), .)

VolzCol0Sig <- sigCol0 %>%
  left_join(., VolzSig) %>%
  mutate(Volz_log2FC = case_when(Volz_log2FC > 0 ~ 'up',
                                      Volz_log2FC < 0 ~ 'down',
                                      is.na(Volz_log2FC) ~ 'no')) %>%
  inner_join(scaleCCol0sig %>% dplyr::select(ID), .)

PauloCol0Sig <- sigCol0 %>%
  left_join(., PauloSig) %>%
  left_join(., PauloHKSig) %>%
  mutate(Paulo_flg22 = case_when(Cluster %in% 1:4 ~ 'up',
                                 Cluster %in% 5:8 ~ 'down',
                                 is.na(Cluster) ~ 'no')) %>%
  mutate(Paulo_bacresp = case_when(logFC > 0 ~ 'up',
                                   logFC < 0 ~ 'down',
                                   is.na(logFC) ~ 'no')) %>%
  dplyr::select(-Cluster, -logFC) %>%
  inner_join(scaleCCol0sig %>% dplyr::select(ID), .)


StringlisFlg22Sig <- sigCol0 %>%
  left_join(sigMatStringlis) %>%
  dplyr::select(ID, Gene, starts_with('flg22')) %>%
  inner_join(scaleCCol0sig %>% dplyr::select(ID), .)

StringlisHKSig <- sigCol0 %>%
  left_join(sigMatStringlis) %>%
  dplyr::select(ID, Gene, starts_with('WCS417')) %>%
  inner_join(scaleCCol0sig %>% dplyr::select(ID), .)

c(sum(scaleCCol0sig$ID == scaleCWERSig$ID),
  sum(scaleCCol0sig$ID == flg22Col0Sig$ID),
  sum(scaleCCol0sig$ID == hkCol0Sig$ID),
  sum(scaleCCol0sig$ID == ironHKLiveSigDay8$ID),
  sum(scaleCCol0sig$ID == ironHKLiveSigDay15$ID),
  ## sum(scaleCCol0sig$ID == ironCol0Sig$ID),
  sum(scaleCCol0sig$ID == CastrilloCol0Sig$ID),
  sum(scaleCCol0sig$ID == VolzCol0Sig$ID),
  sum(scaleCCol0sig$ID == PauloCol0Sig$ID),
  sum(scaleCCol0sig$ID == StringlisFlg22Sig$ID),
  sum(scaleCCol0sig$ID == StringlisHKSig$ID),
  nrow(scaleCCol0sig))

ht_list <- Heatmap(matrix = scaleCCol0sig %>% dplyr::select(contains('_')),
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
  Heatmap(matrix = scaleCWERSig %>% dplyr::select(contains('_')),
          name = 'Scaled Counts',
          column_order = 1 : 12,
          column_split = rep(c('Mock', 'Mock+flg22', 'Non-sup+flg22', 'Sup+flg22'), each = 3),
          show_column_names = FALSE,
          col = colorRampPalette(rev(brewer.pal(n = 10, name = 'Spectral'))[c(-3, -4, -6, -7)])(100),
          column_title_gp = gpar(fontsize = 6),
          use_raster = FALSE) +
  Heatmap(flg22Col0Sig %>% dplyr::select(-ID, -Gene),
          col = c('down' = 'blue', 'no' = 'white', 'up' = 'red'),
          column_names_gp = gpar(fontsize = 5),
          heatmap_legend_param = list(title = 'Col0flg22'),
          cluster_columns = FALSE,
          use_raster = FALSE) +
  Heatmap(PauloCol0Sig %>% dplyr::select(-ID, -Gene, -Paulo_bacresp),
          col = c('down' = 'blue', 'no' = 'white', 'up' = 'red'),
          column_names_gp = gpar(fontsize = 5),
          heatmap_legend_param = list(title = 'Paulo'),
          cluster_columns = FALSE,
          use_raster = FALSE) +
  Heatmap(StringlisFlg22Sig %>% dplyr::select(-ID, -Gene),
          col = c('down' = 'blue', 'no' = 'white', 'up' = 'red'),
          column_names_gp = gpar(fontsize = 5),
          heatmap_legend_param = list(title = 'Stringlis'),
          cluster_columns = FALSE,
          use_raster = FALSE) +
  Heatmap(CastrilloCol0Sig %>% dplyr::select(-ID, -Gene),
          col = c('down' = 'blue', 'no' = 'white', 'up' = 'red'),
          column_names_gp = gpar(fontsize = 5),
          heatmap_legend_param = list(title = 'Castrillo'),
          cluster_columns = FALSE,
          use_raster = FALSE) +
  Heatmap(VolzCol0Sig %>% dplyr::select(-ID, -Gene),
          col = c('down' = 'blue', 'no' = 'white', 'up' = 'red'),
          column_names_gp = gpar(fontsize = 5),
          heatmap_legend_param = list(title = 'Volz'),
          cluster_columns = FALSE,
          use_raster = FALSE) +
  ## Heatmap(ironCol0Sig %>% select(ironbac),
  ##         col = c('bacup' = 'purple', 'bacno' = 'white', 'bacdown' = 'green3'),
  ##         column_names_gp = gpar(fontsize = 5),
  ##         heatmap_legend_param = list(title = 'IronBac'),
  ##         cluster_columns = FALSE,
  ##         use_raster = FALSE) +
  Heatmap(hkCol0Sig %>% dplyr::select(-ID, -Gene),
          col = c('up' = 'red', 'no' = 'white', 'down' = 'blue'),
          column_names_gp = gpar(fontsize = 5),
          heatmap_legend_param = list(title = 'Col0HKlive'),
          cluster_columns = FALSE,
          use_raster = FALSE) +
  Heatmap(PauloCol0Sig %>% dplyr::select(-ID, -Gene, -Paulo_flg22),
          col = c('down' = 'blue', 'no' = 'white', 'up' = 'red'),
          column_names_gp = gpar(fontsize = 5),
          heatmap_legend_param = list(title = 'PauloHKlive'),
          cluster_columns = FALSE,
          use_raster = FALSE) +
  ## Heatmap(ironHKLiveSigDay8 %>% dplyr::select(-3:-6, -ID, -Gene),
  ##         col = c('bacup' = 'red', 'bacno' = 'white', 'bacdown' = 'blue'),
  ##         column_names_gp = gpar(fontsize = 5),
  ##         heatmap_legend_param = list(title = 'IronBac'),
  ##         cluster_columns = FALSE,
  ##         use_raster = FALSE) +
  Heatmap(ironHKLiveSigDay15 %>% dplyr::select(-3:-6,-8:-10, -ID, -Gene),
          col = c('bacup' = 'red', 'bacno' = 'white', 'bacdown' = 'blue'),
          column_names_gp = gpar(fontsize = 5),
          heatmap_legend_param = list(title = 'IronBac'),
          cluster_columns = FALSE,
          use_raster = FALSE) +
Heatmap(StringlisHKSig %>% dplyr::select(-ID, -Gene),
        col = c('bacup' = 'red', 'bacno' = 'white', 'bacdown' = 'blue'),
        column_names_gp = gpar(fontsize = 5),
        heatmap_legend_param = list(title = 'IronBac'),
        cluster_columns = FALSE,
        use_raster = FALSE)
  ## Heatmap(ironCol0Sig %>% select(ironresp),
  ##         col = c('respup' = 'purple', 'respno' = 'white', 'respdown' = 'green3'),
  ##         column_names_gp = gpar(fontsize = 5),
  ##         heatmap_legend_param = list(title = 'IronResp'),
  ##         cluster_columns = FALSE,
  ##         use_raster = FALSE)

## filePrefix <- 'kmeans10_heatmap_WER_Col02'
## filePrefix <- 'kmeans10_heatmap_WER_Col02_Iron2'
## filePrefix <- 'kmeans10_heatmap_WER_Col02_flg22'
filePrefix <- 'kmeans10_heatmap_WER_Col02_flg22_IronDay15'

pdf(paste0(filePrefix, '.pdf'))
draw(ht_list)
dev.off()

system(paste0('convert -density 1200 ', paste0(filePrefix, '.pdf'), ' ', paste0(filePrefix, '.jpg')))
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~HK vs live~~~~~~~~~~~~~~~~~~~~~~~
library('eulerr')
library('VennDiagram')

kmeansRes <- read_csv('kmeans10_1stadd.csv') %>%
  transmute_all(list(~ if_else(is.na(.), FALSE, .)))

fcsig <- kmeansRes %>%
  dplyr::select(ends_with('FoldChange')) %>%
  transmute_all(list(~ case_when(. > log2(1.5) ~ 1,
                                 . < -log2(1.5) ~ -1,
                                 TRUE ~ 0)))

padjsig <- kmeansRes %>%
  dplyr::select(ends_with('padj')) %>%
  `<`(0.05) %>%
  as_tibble

hkCol0SigVenn <- (padjsig * fcsig) %>%
  as_tibble %>%
  setNames(names(.) %>% substr(., start = 1, stop = nchar(.) - 5)) %>%
  transmute_at(.var = vars(contains('vs')),
               list(~ case_when(. == -1 ~ 'bacdown',
                                . == 0 ~'bacno',
                                . == 1 ~ 'bacup'))) %>%
  dplyr::rename(Nonsupp = SynCom33_vs_HKSynCom33,
                Supp = SynCom35_vs_HKSynCom35) %>%
  dplyr::select(Nonsupp, Supp) %>%
  mutate_all(list(~str_detect(., 'down|up'))) %>%
  dplyr::mutate(ID = kmeansRes$ID,
                cl = kmeansRes$cl)

ironHKLiveVennDay8 <- sigMatIronDay8 %>%
  ## dplyr::select(-1:-4, -6:-8, -ID) %>%
  dplyr::select(-1:-4, -ID) %>%
  transmute(IronDay8 = apply(., 1, function(x) {str_detect(x, 'bacdown|bacup') %>% any})) %>%
  mutate(ID = sigMatIronDay8$ID)

ironHKLiveVennDay15 <- sigMatIronDay15 %>%
  ## dplyr::select(-1:-4, -6:-8, -ID) %>%
  dplyr::select(-1:-4, -ID) %>%
  transmute(IronDay15 = apply(., 1, function(x) {str_detect(x, 'bacdown|bacup') %>% any})) %>%
  mutate(ID = sigMatIronDay15$ID)

PauloHKLiveVenn <- PauloHKSig %>%
  dplyr::mutate(Paulo = TRUE) %>%
  dplyr::select(-logFC)

mergeHKVenn <- hkCol0SigVenn %>%
  full_join(ironHKLiveVennDay8) %>%
  full_join(ironHKLiveVennDay15) %>%
  mutate(Gene = ID %>% strsplit(split = '.', fixed = TRUE) %>% sapply('[[', 1)) %>%
  full_join(PauloHKLiveVenn) %>%
  dplyr::filter(!is.na(cl)) %>%
  mutate_all(~ifelse(is.na(.), FALSE, .))

## output venn
allVenn <- foreach (i = 1:10, .combine = inner_join) %do% {

  eachVenn <- mergeHKVenn %>%
    filter(cl == i) %>%
    dplyr::select(IronDay15, Nonsupp, Supp, Paulo) %>%
    euler %>%
    .$original.values %>%
    as.data.frame %>%
    rownames_to_column('ID') %>%
    magrittr::set_colnames(c('ID', paste0('cluster', i)))
}

allVenn %<>%
  mutate(., clusterAll = allVenn %>% dplyr::select(-ID) %>% rowSums)
write_csv(allVenn, 'hk_living_veen.csv')

## plot
foreach (i = 1:10) %do% {

  pdf(paste0('iron_onlysuf_day8_venn/cluster', i, '.pdf'))

  ## mergeHKVenn %>%
  ##   filter(cl == i) %>%
  ##   dplyr::select(Iron, Nonsupp, Supp, Paulo) %>%
  ##   venneuler %>%
  ##   plot

  mergeHKVenn %>%
    ## filter(cl == i) %>%
    dplyr::select(IronDay8, Nonsupp, Supp, Paulo) %>%
    euler(shape = 'ellipse') %>%
    plot(quantities = TRUE,
         labels = list(font = 4),
         fill = c('#ffff33', '#F1696D', '#21BDC3', '#7CAE00'))

  ## mergeHKVenn %>%
  ##   filter(cl == i) %>%
  ##   dplyr::select(IronDay15, Nonsupp, Supp, Paulo) %>% {
  ##     vennList <- list(IronDay15 = which(.$IronDay15),
  ##                      Nonsupp = which(.$Nonsupp),
  ##                      Supp = which(.$Supp),
  ##                      Paulo = which(.$Paulo))
  ##     return(vennList)} %>%
  ##   venn.diagram(filename = paste0('iron_venn3/cluster', i, '.pdf'),
  ##                lwd = 2,
  ##                lty = 'blank',
  ##                fill = c('#ffff33', '#F1696D', '#21BDC3', '#7CAE00'),
  ##                cex = .6,
  ##                fontface = 'bold',
  ##                fontfamily = 'sans',
  ##                cat.cex = 0.6,
  ##                cat.fontface = 'bold',
  ##                cat.default.pos = 'outer',
  ##                cat.fontfamily = 'sans')

  dev.off()
}

## common GO
comHKVeen <- list(
  Nonsupp = mergeHKVenn %>%
    dplyr::filter(Nonsupp) %>%
    .$Gene,
  Supp = mergeHKVenn %>%
    dplyr::filter(Supp) %>%
    .$Gene,
  Paulo = mergeHKVenn %>%
    dplyr::filter(Paulo) %>%
    .$Gene,
  IronDay8 = mergeHKVenn %>%
    dplyr::filter(IronDay8) %>%
    .$Gene,
  ## IronDay15 = mergeHKVenn %>%
  ##   dplyr::filter(IronDay15) %>%
  ##   .$Gene,
  Nonsupp_Supp = mergeHKVenn %>%
    dplyr::filter(Nonsupp, Supp) %>%
    .$Gene,
  Nonsupp_Supp_Paulo = mergeHKVenn %>%
    dplyr::filter(Nonsupp, Supp, Paulo) %>%
    .$Gene,
  ## Nonsupp_Supp_Paulo_IronDay8_IronDay15 = mergeHKVenn %>%
  ##   dplyr::filter(Nonsupp, Supp, Paulo, IronDay8, IronDay15) %>%
  ##   .$Gene
  ## Nonsupp_Supp_Paulo_IronDay15 = mergeHKVenn %>%
  ##   dplyr::filter(Nonsupp, Supp, Paulo, IronDay15) %>%
  ##   .$Gene
  Nonsupp_Supp_Paulo_IronDay8 = mergeHKVenn %>%
    dplyr::filter(Nonsupp, Supp, Paulo, IronDay8) %>%
    .$Gene
) %>%
  compareCluster(geneCluster = .,
                 fun = 'enrichGO',
                 OrgDb = 'org.At.tair.db',
                 keyType= 'TAIR',
                 ont = 'BP',
                 universe = keys(org.At.tair.db),
                 pAdjustMethod = 'BH',
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 1
                 ## pvalueCutoff=0.05,
                 ## qvalueCutoff=0.1
                 )

dotplot(comHKVeen, showCategory = 30, font.size = 8)
ggsave('common_HKvsLiving.pdf', height = 20)
write_csv(as.data.frame(comHKVeen), 'common_HKvsLiving.csv')

##~~~~~~~~~~~~~~~~~~~~~~~~~~core response genes~~~~~~~~~~~~~~~~~~~~~~~
library('GO.db')

coreHKLive <- mergeHKVenn %>%
  dplyr::filter(Nonsupp, Supp, Paulo, IronDay8) %>%
  .$Gene

xx <- as.list(org.At.tairGO[coreHKLive])

CollpaseEachGO <- function(eachGOList, geneID) {

  require('foreach')
  require('tidyverse')

  res <- foreach(i = seq_along(eachGOList), .combine = rbind) %do% {
    return(unlist(eachGOList[[i]]))
  } %>%
    as_tibble %>%
  mutate(GeneID = geneID)

  return(res)
}

allGOAnno <- Term(GOTERM) %>%
  unlist %>%
  as.data.frame %>%
  rownames_to_column('GOID') %>%
  set_colnames(c('GOID', 'Term')) %>%
  as_tibble

coreHKLiveMat <- foreach(i = seq_along(xx), .combine = bind_rows) %do% {
  res <- CollpaseEachGO(xx[[i]], names(xx)[i])
} %>%
  dplyr::filter(Ontology %in% c('BP')) %>%
  dplyr::select(-value, -Evidence, -Ontology) %>%
  dplyr::distinct(GOID, GeneID) %>%
  dplyr::inner_join(allGOAnno) %>%
  dplyr::group_by(GeneID) %>%
  dplyr::summarize_all(~paste(., collapse = ';'))

read_csv('kmeans10_1stadd_sig.csv') %>%
  mutate_all(~ifelse(is.na(.), '', .)) %>%
  dplyr::select(ID, Gene, Description, cl) %>%
  dplyr::mutate(GeneID = ID %>% strsplit('.', fixed = TRUE) %>% sapply('[[', 1)) %>%
  dplyr::inner_join(coreHKLiveMat) %>%
  write_csv('hk_living_core.csv')

## compare geneID
comHKVeen %>%
  as.data.frame %>%
  as_tibble %>%
  dplyr::filter(Cluster %in% c('Nonsupp_Supp_Paulo_IronDay8')) %>%
  .$geneID %>%
  strsplit(split = '/', fixed = TRUE) %>%
  unlist %>%
  unique %>%
  length
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## select general response plot
comVeenSelect <- comHKVeen
comVeenSelect@compareClusterResult <-
  c('GO:0010054', 'GO:0010015', 'GO:0098754',
    'GO:0001666', 'GO:0015698', 'GO:0071453',
    'GO:0010200', 'GO:0050832', 'GO:0009862',
    'GO:0009991') %>%
  {comVeenSelect@compareClusterResult$ID %in% .} %>%
  which %>%
  comVeenSelect@compareClusterResult[., ]

dotplot(comVeenSelect, showCategory = 30)
ggsave('common_HKvsLiving_selectGO.pdf', width = 20)

## general defense
defenseTerms <- c('GO:0010243', 'GO:0010200', 'GO:0009697',
                  'GO:0009696', 'GO:0050832', 'GO:0009863',
                  'GO:0071446', 'GO:0009626', 'GO:0034050',
                  'GO:0009627', 'GO:0002679', 'GO:0045730',
                  'GO:0009753', 'GO:0009723', 'GO:0060548',
                  'GO:0000302', 'GO:0043069', 'GO:0036294',
                  'GO:0010112', 'GO:0010337', 'GO:0012501',
                  'GO:0009862')

defenseGenes <- comVeen %>%
  as.data.frame %>%
  as_tibble %>%
  dplyr::filter(ID %in% defenseTerms) %>%
  dplyr::group_by(Cluster) %>%
  dplyr::summarise(geneID = paste(geneID, collapse = '/')) %>% {
    uniGene <- .$geneID %>%
      strsplit(split = '/', fixed = TRUE) %>%
      sapply(function(x) {x %>% unique %>% paste(collapse = '/')})

    mutate(., geneID = uniGene)
  } %>% {
    len <- .$geneID %>%
      strsplit(split = '/', fixed = TRUE) %>%
      sapply(length)

    mutate(., No = len)
  }

write_csv(defenseGenes, 'Nonsupp_Supp_Paulo_Iron_defenseGenes.csv')

## iron test
sigMatIronFull <- full_join(sigMatIronDay8 %>%
          dplyr::select(-1:-4) %>%
          transmute_at(.var = vars(contains('vs')),
                       list(~ case_when(. %in% c('bacup', 'bacdown') ~ TRUE,
                                        . %in% 'bacno' ~ FALSE))) %>%
          mutate(ID = sigMatIronDay8$ID)
         ,
          sigMatIronDay15 %>%
          dplyr::select(-1:-4) %>%
          transmute_at(.var = vars(contains('vs')),
                       list(~ case_when(. %in% c('bacup', 'bacdown') ~ TRUE,
                                        . %in% 'bacno' ~ FALSE))) %>%
          mutate(ID = sigMatIronDay15$ID)) %>%
  mutate_all(~ifelse(is.na(.), FALSE, .)) %>%
  dplyr::rename(Col0_Suf_Day8 = Col0_FeEDTA_Live_vs_Col0_FeEDTA_HK_Day8,
                Col0_Def_Day8 = Col0_FeCl3_Live_vs_Col0_FeCl3_HK_Day8,
                f6h1_Suf_Day8 = f6h1_FeEDTA_Live_vs_f6h1_FeEDTA_HK_Day8,
                f6h1_Def_Day8 = f6h1_FeCl3_Live_vs_f6h1_FeCl3_HK_Day8,
                Col0_Suf_Day15 = Col0_FeEDTA_Live_vs_Col0_FeEDTA_HK_Day15,
                Col0_Def_Day15 = Col0_FeCl3_Live_vs_Col0_FeCl3_HK_Day15,
                f6h1_Suf_Day15 = f6h1_FeEDTA_Live_vs_f6h1_FeEDTA_HK_Day15,
                f6h1_Def_Day15 = f6h1_FeCl3_Live_vs_f6h1_FeCl3_HK_Day15)

pdf('iron_only_venn/iron_day8.pdf')
sigMatIronFull %>%
  dplyr::select(-6:-9, -ID) %>%
  euler(shape = 'ellipse') %>%
  plot(quantities = TRUE,
       labels = list(font = 4))
dev.off()
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~flg22 Venn~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library('eulerr')
library('magrittr')

kmeansRes <- read_csv('kmeans10_1stadd.csv') %>%
  mutate_all(~ifelse(is.na(.), '', .))

fcsig <- kmeansRes %>%
  dplyr::select(ends_with('FoldChange')) %>%
  transmute_all(list(~ case_when(. > log2(1.5) ~ 1,
                                 . < -log2(1.5) ~ -1,
                                 TRUE ~ 0)))

padjsig <- kmeansRes %>%
  dplyr::select(ends_with('padj')) %>%
  `<`(0.05) %>%
  as_tibble %>%
  transmute_all(list(~ if_else(is.na(.), FALSE, .)))

flg22Col0SigVenn <- (padjsig * fcsig) %>%
  as_tibble %>%
  setNames(names(.) %>% substr(., start = 1, stop = nchar(.) - 5)) %>%
  transmute_at(.var = vars(contains('vs')),
               list(~ case_when(. == -1 ~ 'bacdown',
                                . == 0 ~'bacno',
                                . == 1 ~ 'bacup'))) %>%
  dplyr::rename(Mock = Mock_Flg22_vs_Mock,
                HKNonsupp = HKSynCom33_Flg22_vs_HKSynCom33,
                HKSupp = HKSynCom35_Flg22_vs_HKSynCom35,
                Nonsupp = SynCom33_Flg22_vs_SynCom33,
                Supp = SynCom35_Flg22_vs_SynCom35) %>%
  dplyr::select(Mock, HKNonsupp, HKSupp, Nonsupp, Supp) %>%
  mutate_all(list(~str_detect(., 'bacdown|bacup'))) %>%
  dplyr::mutate(ID = kmeansRes$ID,
                cl = kmeansRes$cl)

flg22PauloSigVenn <- PauloSig %>%
  mutate(Paulo = case_when(Cluster %in% 1:4 ~ 'up',
                           Cluster %in% 5:8 ~ 'down',
                           is.na(Cluster) ~ 'no')) %>%
  mutate(Paulo = Paulo %>% str_detect('down|up')) %>%
  dplyr::select(Gene, Paulo)

flg22StringlisSigVenn <- sigMatStringlis %>%
  dplyr::rename(flg22Pa_0p5h = flg22Pa_0p5h_vs_Mock_0p5h,
                flg22Pa_1h = flg22Pa_1h_vs_Mock_1h,
                flg22Pa_3h = flg22Pa_3h_vs_Mock_3h,
                flg22Pa_6h = flg22Pa_6h_vs_Mock_6h,
                flg22417_0p5h = flg22417_0p5h_vs_Mock_0p5h,
                flg22417_1h = flg22417_1h_vs_Mock_1h,
                flg22417_3h = flg22417_3h_vs_Mock_3h,
                flg22417_6h = flg22417_6h_vs_Mock_6h) %>%
  dplyr::select(1:8) %>%
  mutate_all(list(~str_detect(., 'down|up'))) %>%
  dplyr::mutate(ID = sigMatStringlis$ID)

mergeflg22Venn <- flg22Col0SigVenn %>%
  full_join(flg22StringlisSigVenn) %>%
  mutate(Gene = ID %>% strsplit(split = '.', fixed = TRUE) %>% sapply('[[', 1)) %>%
  full_join(flg22PauloSigVenn) %>%
  mutate(Mock_HK = apply(.[, 1:4], 1, any)) %>%
  dplyr::filter(!is.na(cl)) %>%
  mutate_all(~ifelse(is.na(.), FALSE, .))

## output venn
allVenn <- foreach (i = 1:10, .combine = inner_join) %do% {

  eachVenn <- mergeflg22Venn %>%
    filter(cl == i) %>%
    dplyr::select(Mock_HK, Nonsupp, Paulo, flg22Pa_6h, flg22417_6h) %>%
    euler %>%
    .$original.values %>%
    as.data.frame %>%
    rownames_to_column('ID') %>%
    magrittr::set_colnames(c('ID', paste0('cluster', i)))
}

allVenn %<>%
  mutate(., clusterAll = allVenn %>% dplyr::select(-ID) %>% rowSums)
write_csv(allVenn, 'flg22_veen.csv')

foreach (i = 1:10) %do% {

  pdf(paste0('flg22_venn/cluster', i, '.pdf'))

  mergeflg22Venn %>%
    filter(cl == i) %>%
    dplyr::select(Mock_HK, Nonsupp, Paulo, flg22Pa_6h, flg22417_6h) %>%
    euler(shape = 'ellipse') %>%
    plot(quantities = TRUE,
         labels = list(font = 2))

  dev.off()
}

## common GO
comflg22Veen <- list(
  Mock_HK = mergeflg22Venn %>%
    dplyr::filter(Mock_HK) %>%
    .$Gene,
  Nonsupp = mergeflg22Venn %>%
    dplyr::filter(Nonsupp) %>%
    .$Gene,
  Paulo = mergeflg22Venn %>%
    dplyr::filter(Paulo) %>%
    .$Gene,
  flg22Pa_6h = mergeflg22Venn %>%
    dplyr::filter(flg22Pa_6h) %>%
    .$Gene,
  flg22417_6h = mergeflg22Venn %>%
    dplyr::filter(flg22417_6h) %>%
    .$Gene,
  ## Nonsupp_Supp = mergeflg22Venn %>%
  ##   dplyr::filter(Nonsupp, Supp) %>%
  ##   .$Gene,
  ## Nonsupp_Supp_Paulo = mergeflg22Venn %>%
  ##   dplyr::filter(Nonsupp, Supp, Paulo) %>%
  ##   .$Gene,
  MockHK_Nonsupp_Paulo_flg22Pa6h_flg224176h = mergeflg22Venn %>%
    dplyr::filter(Mock_HK, Nonsupp, Paulo, flg22Pa_6h, flg22417_6h) %>%
    .$Gene
) %>%
  compareCluster(geneCluster = .,
                 fun = 'enrichGO',
                 OrgDb = 'org.At.tair.db',
                 keyType= 'TAIR',
                 ont = 'BP',
                 universe = keys(org.At.tair.db),
                 pAdjustMethod = 'BH',
                 pvalueCutoff=0.05,
                 qvalueCutoff=0.1)

dotplot(comVeen, showCategory = 20, font.size = 8)
ggsave('common_flg22.pdf', height = 20)
write_csv(as.data.frame(comVeen), 'common_flg22.csv')

pdf('flg22_response_DEGs.pdf')
flg22Col %>% transmute_at(.var = vars(-ID),
                          list(~ case_when(str_detect(., 'down|up') ~ TRUE,
                                           str_detect(., 'no') ~ FALSE))) %>%
  set_colnames(c('Mock+flg22 vs. Mock', 'HK_NS3+flg22 vs. HK_NS3', 'HK_S3+flg22 vs. HK_S3', 'NS3+flg22 vs. NS3', 'S3+flg22 vs. S3')) %>%
  euler %>%
  plot(quantities = TRUE,
       labels = list(font = 4),
       fill = c(NA, '#d9d9d9', '#d9d9d9', '#F1696D'))
dev.off()
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#######################################################################

##########################Col0 agar soil##############################
library('tidyverse')
library('ComplexHeatmap')
library('limma')
library('DESeq2')
library('RColorBrewer')
library('clusterProfiler')

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


###########################compare Paulo##############################
library('tidyverse')
library('magrittr')
library('foreach')
library('corrplot')
library('RColorBrewer')
library('org.At.tair.db')
library('clusterProfiler')

basepath <- '/extDisk1/RESEARCH'

basepath %>%
  file.path('MPIPZ_KaWai_RNASeq/results/removeZero') %>%
  setwd

kmeansResSig <- read_csv('kmeans10_1stadd_sig.csv') %>%
  dplyr::select(ID, cl)

Paulo <- read_delim('Paulo_RNASeq_flg22.txt', delim = '\t')

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~Jaccard similarity~~~~~~~~~~~~~~~~~~~~~~~
Col0PauloJS <- foreach (i = 1:10, .combine = rbind) %do% {

  Col0Each <- kmeansResSig %>%
    filter(cl == i) %>%
    .$ID %>%
    substr(start = 1, stop = nchar(.) - 2) %>%
    unique

  JacSim <- foreach(j = 1:8, .combine = c)  %do% {
    PauloEach <- Paulo %>%
      filter(Cluster == j) %>%
      .$Gene %>%
      unique

    eachJacSim <- length(intersect(Col0Each, PauloEach)) / length(union(Col0Each, PauloEach))

    return(eachJacSim)

  } %>%
  set_names(paste0('Paulo', 1:8))

  return(JacSim)
} %>%
  set_rownames(paste0('Ka-Wai', 1:10)) %>%
  round(digits = 3)

Col0PauloJS <- apply(Col0PauloJS, 1:2, function(x){ifelse(x < 0.0001, 0, x)})

pdf('KaWai_Paulo_JacSim.pdf')
cols <- rev(brewer.pal(n = 11, name = 'RdBu'))[c(-3, -4, -5, -7, -8, -9)]
colorRampPalette(cols)(50) %>%
   corrplot(Col0PauloJS, col = ., cl.lim = c(0, 0.5))
dev.off()
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~cluster~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Col0List <- foreach (i = 1:10) %do% {

  kmeansResSig %>%
    filter(cl == i) %>%
    .$ID %>%
    substr(start = 1, stop = nchar(.) - 2) %>%
    unique

} %>%
  set_names(paste0('Ka-Wai', 1:10))

PauloList <- split(Paulo$Gene, Paulo$Cluster) %>%
  lapply(unique) %>%
  set_names(paste0('Paulo', 1:8))

Col0PauloList <- list(
  `Ka-Wai1-Paulo1` = c(1, 1),
  `Ka-Wai1-Paulo2` = c(1, 2),
  `Ka-Wai1-Paulo4` = c(1, 4),
  `Ka-Wai2-Paulo5` = c(2, 5),
  `Ka-Wai2-Paulo8` = c(2, 8),
  `Ka-Wai3-Paulo1` = c(3, 1),
  `Ka-Wai5-Paulo6` = c(5, 6),
  `Ka-Wai5-Paulo8` = c(5, 8)) %>%
  lapply(function(x) {
    intersect(Col0List[[x[1]]], PauloList[[x[2]]])
  })

resList <- c(Col0List[1],
             PauloList[c(1, 2, 4)],
             Col0PauloList[1:3],
             Col0List[2],
             PauloList[c(5, 8)],
             Col0PauloList[4:5],
             Col0List[3],
             PauloList[1],
             Col0PauloList[6],
             Col0List[5],
             PauloList[c(6, 8)],
             Col0PauloList[7:8])

1:7
8:12
13:15
16:20
goBP <- compareCluster(geneCluster = resList[13:15],
                       fun = 'enrichGO',
                       OrgDb = 'org.At.tair.db',
                       keyType= 'TAIR',
                       ont = 'BP',
                       universe = keys(org.At.tair.db),
                       pAdjustMethod = 'BH',
                       pvalueCutoff=0.01,
                       qvalueCutoff=0.1)

dotplot(goBP, showCategory = 20)
ggsave('compare_Ka-Wai_Paulo_4.pdf', width = 15, height = 12)

write_csv(as.data.frame(goBP), 'compare_Ka-Wai3_Paulo1.csv')

## Col0 Paulo common defense response genes
goBPdefense <- compareCluster(geneCluster = resList[13:15],
                              fun = 'enrichGO',
                              OrgDb = 'org.At.tair.db',
                              keyType= 'TAIR',
                              ont = 'BP',
                              universe = keys(org.At.tair.db),
                              pAdjustMethod = 'BH',
                              pvalueCutoff=0.01,
                              qvalueCutoff=0.1)
defenseTerms <- c('GO:0002679', 'GO:0045730', 'GO:0010200',
                  'GO:0010243', 'GO:0009753', 'GO:0050832',
                  'GO:0009723', 'GO:0009867', 'GO:0071395',
                  'GO:0002252', 'GO:0009751', 'GO:0012501',
                  'GO:0010941', 'GO:0010363', 'GO:0008219',
                  'GO:0012501', 'GO:0009626', 'GO:0034050')

defenseGenes <- goBPdefense %>%
  as.data.frame %>%
  as_tibble %>%
  dplyr::filter(ID %in% defenseTerms) %>%
  dplyr::group_by(Cluster) %>%
  dplyr::summarise(geneID = paste(geneID, collapse = '/')) %>% {
    uniGene <- .$geneID %>%
      strsplit(split = '/', fixed = TRUE) %>%
      sapply(function(x) {x %>% unique %>% paste(collapse = '/')})

    mutate(., geneID = uniGene)
  } %>% {
    len <- .$geneID %>%
      strsplit(split = '/', fixed = TRUE) %>%
      sapply(length)

    mutate(., No = len)
  }

write_csv(defenseGenes, 'Col0_Paulo_defenseGenes.csv')

## compare defense responsive genes
commDef <- read_csv('Col0_Paulo_defenseGenes.csv')
genDef <- read_csv('Nonsupp_Supp_Paulo_Iron_defenseGenes.csv')

lapply(1:8, function(x) {
  intersect(commDef$geneID[3] %>%
            strsplit(split = '/', fixed = TRUE) %>%
            unlist,
            genDef$geneID[x] %>%
            strsplit(split = '/', fixed = TRUE) %>%
            unlist)
}) %>%
  set_names(genDef$Cluster)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
######################################################################
