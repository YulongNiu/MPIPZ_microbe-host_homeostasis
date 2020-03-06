###########################replot heatmap##############################
library('tidyverse')
library('ComplexHeatmap')
library('limma')
library('DESeq2')
library('RColorBrewer')

setwd('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/removeZero')

load('degres_condi_Mock.RData')

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~select DEGs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
wholeDEG <- read_csv('eachGroup_vs_Mock_k.csv')
kmeansRes <- read_csv('kmeans10.csv') %>%
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
  write_csv('kmeans10_sig.csv')
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

colnames(sigMat) <- c('Mock+flg22 vs. Mock',
                      'SynCom33+flg22 vs. Mock',
                      'SynCom35+flg22 vs. Mock',
                      'SynCom33+flg22 vs Mock+flg22',
                      'SynCom35+flg22 vs. Mock+flg22',
                      'SynCom33+flg22 vs. SynCom35+flg22')

ht_list <- Heatmap(matrix = scaleC %>% select(contains('_')),
                   name = 'Scaled Counts',
                   ## row_order = order(scaleC$cl) %>% rev,
                   row_split = scaleC$cl,
                   row_gap = unit(2, "mm"),
                   column_order = 1 : 12,
                   column_split = rep(c('Mock', 'Mock+flg22', 'Non-sup+flg22', 'Sup+flg22'), each = 3),
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

filePrefix <- 'kmeans10_heatmap_sig_DEG2'

pdf(paste0(filePrefix, '.pdf'))
draw(ht_list)
dev.off()

system(paste0('convert -density 1200 ', paste0(filePrefix, '.pdf'), ' ', paste0(filePrefix, '.jpg')))
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~box plot~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
meanFlg22 <- function(v) {

  require('magrittr')

  res <- v %>%
    split(rep(1 : 4, each = 3)) %>%
    sapply(mean, na.rm = TRUE)

  return(res)
}

## kmeansRes <- read_csv('kmeans10.csv') %>%
##   select(ID, cl)

kmeansResSig <- read_csv('kmeans10_sig.csv') %>%
  select(ID, cl)

sampleN <- c('Mock', 'Mock_flg22', 'Nonsupp_flg22', 'Supp_flg22')

boxplotData <- rldData %>%
  t %>%
  scale %>%
  t %>%
  apply(1, meanFlg22) %>%
  t %>%
  set_colnames(sampleN) %>%
  as.data.frame %>%
  rownames_to_column('ID') %>%
  as_tibble %>%
  inner_join(kmeansResSig)

for (i in 1:10) {
  boxplotData %>%
    filter(cl == i) %>%
    select(-ID, -cl) %>%
    gather(key = 'Conditions', value = 'ScaleCounts') %>%
    mutate(Group = case_when(
             str_detect(Conditions, '^Mock|^HK') ~ 'Mock',
             str_detect(Conditions, '^Nonsupp') ~ 'Nonsupp',
             str_detect(Conditions, '^Supp') ~ 'Supp',
           )) %>%
    mutate(Conditions = Conditions %>% factor(levels = sampleN)) %>%
    mutate(Group = Group %>% factor(levels = c('Mock', 'Nonsupp', 'Supp'))) %>%
    ggplot(aes(x = Group, y = ScaleCounts, fill = Conditions)) +
    geom_boxplot(position = position_dodge2(preserve = 'single')) +
    ## scale_fill_manual(values = c(rep(NA, 6), rep('#377eb8', 2), rep('#e41a1c', 2))) +
    scale_fill_manual(values = c(rep(NA, 2), '#377eb8', '#e41a1c')) +
    ylim(-2, 2) +
    ylab('Scaled counts') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90),
          plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
          legend.text.align = 0,
          axis.text = element_text(size = 13),
          axis.title = element_text(size = 14),
          legend.text=element_text(size= 13),
          legend.title = element_text(size = 14))

  ggsave(paste0('boxplot/kmeans10_boxplot', i, '.pdf'))
  ggsave(paste0('boxplot/kmeans10_boxplot', i, '.jpeg'))
}
#######################################################################


