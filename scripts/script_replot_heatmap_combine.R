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
  file.path('MPIPZ_CJ_RNASeq/results/eachGroup_mergeDay8_deg_sig.csv') %>%
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

sigMatIron <- (padjsigIron * fcsigIron) %>%
  as_tibble %>%
  setNames(names(.) %>% substr(., start = 1, stop = nchar(.) - 5)) %>%
  transmute_at(.var = vars(contains('vs')),
               list(~ case_when(. == -1 ~ 'bacdown',
                                . == 0 ~'bacno',
                                . == 1 ~ 'bacup'))) %>%
  mutate(ID = ironHKLive$ID)

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

ironHKLiveSig <- sigCol0 %>%
  left_join(., sigMatIron) %>%
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

c(sum(scaleCCol0sig$ID == scaleCWERSig$ID),
  sum(scaleCCol0sig$ID == flg22Col0Sig$ID),
  sum(scaleCCol0sig$ID == hkCol0Sig$ID),
  sum(scaleCCol0sig$ID == ironHKLiveSig$ID),
  ## sum(scaleCCol0sig$ID == ironCol0Sig$ID),
  sum(scaleCCol0sig$ID == CastrilloCol0Sig$ID),
  sum(scaleCCol0sig$ID == VolzCol0Sig$ID),
  sum(scaleCCol0sig$ID == PauloCol0Sig$ID),
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
  Heatmap(ironHKLiveSig %>% dplyr::select(-2:-5, -ID, -Gene),
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
filePrefix <- 'kmeans10_heatmap_WER_Col02_flg22_Iron'

pdf(paste0(filePrefix, '.pdf'))
draw(ht_list)
dev.off()

system(paste0('convert -density 1200 ', paste0(filePrefix, '.pdf'), ' ', paste0(filePrefix, '.jpg')))
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~HK vs live~~~~~~~~~~~~~~~~~~~~~~~
library('eulerr')
library('VennDiagram')

hkCol0SigVenn <- hkCol0Sig %>%
  mutate_at(vars(-ID, -Gene), list(~str_detect(., 'down|up'))) %>%
  dplyr::rename(Nonsupp = `SynCom33 vs. HKSynCom33`,
                Supp = `SynCom35 vs. HKSynCom35`)

ironHKLiveVenn <- ironHKLiveSig %>%
  dplyr::select(-1:-6) %>%
  transmute(Iron = apply(., 1, function(x) {str_detect(x, 'bacdown|bacup') %>% any})) %>%
  mutate(ID = ironHKLiveSig$ID)

PauloHKLiveVenn <- PauloCol0Sig %>%
  mutate(Paulo = str_detect(Paulo_bacresp, 'down|up')) %>%
  dplyr::select(ID, Paulo)

mergeVenn <- inner_join(hkCol0SigVenn, ironHKLiveVenn) %>%
  inner_join(PauloHKLiveVenn) %>%
  inner_join(scaleCCol0sig %>% dplyr::select(ID, cl))

## output venn
allVenn <- foreach (i = 1:10, .combine = inner_join) %do% {

  eachVenn <- mergeVenn %>%
    filter(cl == i) %>%
    dplyr::select(Iron, Nonsupp, Supp, Paulo) %>%
    euler %>%
    .$original.values %>%
    as.data.frame %>%
    rownames_to_column('ID') %>%
    magrittr::set_colnames(c('ID', paste0('cluster', i)))
}

allVenn %<>%
  mutate(clusterAll = allVenn %>% dplyr::select(-ID) %>% rowSums)
write_csv(allVenn, 'hk_living_veen.csv')

## plot
foreach (i = 1:10) %do% {

  ## pdf(paste0('iron_venn/cluster', i, '.pdf'))

  ## mergeVenn %>%
  ##   filter(cl == i) %>%
  ##   dplyr::select(Iron, Nonsupp, Supp, Paulo) %>%
  ##   venneuler %>%
  ##   plot

  ## mergeVenn %>%
  ##   filter(cl == i) %>%
  ##   dplyr::select(Iron, Nonsupp, Supp, Paulo) %>%
  ##   euler %>%
  ##   plot(quantities = TRUE,
  ##        labels = list(font = 4),
  ##        fill = c('#ffff33', '#F1696D', '#21BDC3', '#7CAE00'))

  mergeVenn %>%
    filter(cl == i) %>%
    dplyr::select(Iron, Nonsupp, Supp, Paulo) %>% {
      vennList <- list(Iron = which(.$Iron),
                       Nonsupp = which(.$Nonsupp),
                       Supp = which(.$Supp),
                       Paulo = which(.$Paulo))
      return(vennList)} %>%
    venn.diagram(filename = paste0('iron_venn2/cluster', i, '.pdf'),
                 lwd = 2,
                 lty = 'blank',
                 fill = c('#ffff33', '#F1696D', '#21BDC3', '#7CAE00'),
                 cex = .6,
                 fontface = 'bold',
                 fontfamily = 'sans',
                 cat.cex = 0.6,
                 cat.fontface = 'bold',
                 cat.default.pos = 'outer',
                 cat.fontfamily = 'sans')

  ## dev.off()
}

## common GO
comVeen <- list(Nonsupp = mergeVenn %>%
                  dplyr::filter(Nonsupp) %>%
                  .$Gene,
                Supp = mergeVenn %>%
                  dplyr::filter(Supp) %>%
                  .$Gene,
                Paulo = mergeVenn %>%
                  dplyr::filter(Paulo) %>%
                  .$Gene,
                Iron = mergeVenn %>%
                  dplyr::filter(Iron) %>%
                  .$Gene,
                Nonsupp_Supp = mergeVenn %>%
                  dplyr::filter(Nonsupp, Supp) %>%
                  .$Gene,
                ## Supp_Paulo = mergeVenn %>%
                ##   dplyr::filter(Supp, Paulo) %>%
                ##   .$Gene,
                ## Supp_Iron = mergeVenn %>%
                ##   dplyr::filter(Supp, Iron) %>%
                ##   .$Gene,
                Nonsupp_Supp_Paulo = mergeVenn %>%
                  dplyr::filter(Nonsupp, Supp, Paulo) %>%
                  .$Gene,
                Nonsupp_Supp_Iron = mergeVenn %>%
                  dplyr::filter(Nonsupp, Supp, Iron) %>%
                  .$Gene,
                Nonsupp_Supp_Paulo_Iron = mergeVenn %>%
                  dplyr::filter(Nonsupp, Supp, Paulo, Iron) %>%
                  .$Gene) %>%
  compareCluster(geneCluster = .,
                 fun = 'enrichGO',
                 OrgDb = 'org.At.tair.db',
                 keyType= 'TAIR',
                 ont = 'BP',
                 universe = keys(org.At.tair.db),
                 pAdjustMethod = 'BH',
                 pvalueCutoff=0.05,
                 qvalueCutoff=0.1)

dotplot(comVeen, showCategory = 40, font.size = 8)
ggsave('common_HKvsLiving.pdf', height = 20)

write_csv(as.data.frame(comVeen), 'common_HKvsLiving.csv')


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
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~flg22 Venn~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library('eulerr')
library('magrittr')


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
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
######################################################################
