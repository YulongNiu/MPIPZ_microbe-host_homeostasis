
## originally by Yulong Niu
## yulong.niu@hotmail.com

#############################output each sigpath########################
setwd('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/geneset/')

library('goseq')
library('GO.db')
library('foreach')
library('doMC')
library('KEGGAPI')
library('BioCycAPI')
library('magrittr')
library('dplyr')
library('tibble')
library('readr')
library('stringr')

registerDoMC(12)

load('athGO.RData')
load('athKEGG.RData')
load('athBioCyc.RData')

kmeansRes <- read_csv('../results/kmeans_10.csv',
                      col_types = cols(Chromosome = col_character()))
savePath <- '/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/pathway/splitpath'

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~choose sig~~~~~~~~~~~~~~~~~~~~~
## padj < 0.05 & |log2FC| > log2(1.5)
fcsig <- kmeansRes %>%
  select(ends_with('FoldChange')) %>%
  transmute_all(list(~ case_when(. > 1 ~ 1,
                                 . < -1 ~ -1,
                                 TRUE ~ 0)))
padjsig <- kmeansRes %>%
  select(ends_with('padj')) %>%
  `<`(0.05) %>%
  as_tibble %>%
  transmute_all(list(~ if_else(is.na(.), FALSE, .)))

sig <- (padjsig * fcsig) %>%
  as_tibble %>%
  mutate(ID = kmeansRes$ID, cl = kmeansRes$cl) %>%
  select(ID, everything()) %>%
  rename(Flg22_vs_Mock = Flg22_vs_Mock_padj,
         Flg22_SynCom33_vs_Mock = Flg22_SynCom33_vs_Mock_padj,
         Flg22_SynCom35_vs_Mock = Flg22_SynCom35_vs_Mock_padj) %>%
  {kanno <- kmeansRes %>%
     dplyr::select(ID : Description, Flg22_vs_Mock_pvalue : Flg22_SynCom35_vs_Mock_log2FoldChange)
     inner_join(kanno, .)} %>%
  mutate(Gene = Gene %>% coalesce('')) %>%
  mutate(Description = Description %>% coalesce(''))
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~select genesets~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
athGO %<>%
  lapply(function(x) {x[x %in% kmeansRes$ID]}) %>%
  {
    l <- sapply(., length) > 0
    .[l]
  }
GOMat <- foreach(i = 1:length(athGO), .combine = rbind) %dopar% {
  eachMat <- cbind(athGO[[i]], names(athGO)[i])
  return(eachMat)
} %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  set_colnames(c('ID', 'Path')) %>%
  as_tibble

GOMat %<>% inner_join(sig)
GOList <- split(GOMat, GOMat$Path)
for(i in seq_along(GOList)) {
  eachpath <- paste0(savePath, '/', names(GOList)[i], '.csv') %>%
    write_csv(GOList[[i]], .)
}

athKEGG %<>%
  lapply(function(x) {x[x %in% kmeansRes$ID]}) %>%
  {
    l <- sapply(., length) > 0
    .[l]
  }
KEGGMat <- foreach(i = 1:length(athKEGG), .combine = rbind) %dopar% {
  eachMat <- cbind(athKEGG[[i]], names(athKEGG)[i])
  return(eachMat)
} %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  set_colnames(c('ID', 'Path')) %>%
  as_tibble

KEGGMat %<>% inner_join(sig)
KEGGList <- split(KEGGMat, KEGGMat$Path)
for(i in seq_along(KEGGList)) {
  eachpath <- paste0(savePath, '/', names(KEGGList)[i], '.csv') %>%
    write_csv(KEGGList[[i]], .)
}


athBioCyc %<>%
  lapply(function(x) {x[x %in% kmeansRes$ID]}) %>%
  {
    l <- sapply(., length) > 0
    .[l]
  }
BioCycMat <- foreach(i = 1:length(athBioCyc), .combine = rbind) %dopar% {
  eachMat <- cbind(athBioCyc[[i]], names(athBioCyc)[i])
  return(eachMat)
} %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  set_colnames(c('ID', 'Path')) %>%
  as_tibble

BioCycMat %<>% inner_join(sig)
BioCycList <- split(BioCycMat, BioCycMat$Path)
for(i in seq_along(BioCycList)) {
  eachpath <- paste0(savePath, '/', names(BioCycList)[i], '.csv') %>%
    write_csv(BioCycList[[i]], .)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#######################################################################
