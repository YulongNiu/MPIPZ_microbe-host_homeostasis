#######################GO analysis############################
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

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~select DEGs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/removeZero')

wholeDEG <- read_csv('eachGroup_vs_Mock_k_1stadd.csv')

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
  setNames(., substring(colnames(.), first = 1, last = nchar(colnames(.)) - 5)) %>%
  bind_cols(wholeDEG %>% select(ID, Length))
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~GO~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
registerDoMC(12)

setwd('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/geneset')

load('athGO.RData')
load('athKEGG.RData')
load('athBioCyc.RData')

savepath <- '/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/removeZero/geneset_1stadd/DEGs'

##~~~~~~~~~~~~~~~select genesets~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
athGO %<>%
  lapply(function(x) {x[x %in% sigMat$ID]}) %>%
  {
    l <- sapply(., length) > 0
    .[l]
  }
GOMat <- foreach(i = 1:length(athGO), .combine = rbind) %dopar% {
  eachMat <- cbind(athGO[[i]], names(athGO)[i])
  return(eachMat)
} %>% as.data.frame

athKEGG %<>%
  lapply(function(x) {x[x %in% sigMat$ID]}) %>%
  {
    l <- sapply(., length) > 0
    .[l]
  }
KEGGMat <- foreach(i = 1:length(athKEGG), .combine = rbind) %dopar% {
  eachMat <- cbind(athKEGG[[i]], names(athKEGG)[i])
  return(eachMat)
} %>% as.data.frame


athBioCyc %<>%
  lapply(function(x) {x[x %in% sigMat$ID]}) %>%
  {
    l <- sapply(., length) > 0
    .[l]
  }
BioCycMat <- foreach(i = 1:length(athBioCyc), .combine = rbind) %dopar% {
  eachMat <- cbind(athBioCyc[[i]], names(athBioCyc)[i])
  return(eachMat)
} %>% as.data.frame
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~whole cluster gene-set~~~~~~~~~~~~~~~~~
cond <- colnames(sigMat) %>%
  {.[-which(. %in% c('ID', 'Length'))]}

for (i in cond) {

  ## up
  degVec <- (sigMat[, i] > 0) %>%
    as.integer %>%
    set_names(sigMat$ID)
  pwfup <- nullp(degVec, bias.data = sigMat$Length)

  ## down
  degVec <- (sigMat[, i] < 0) %>%
    as.integer %>%
    set_names(sigMat$ID)
  pwfdown <- nullp(degVec, bias.data = sigMat$Length)

  ## GO
  GOTestWithCatUp <- goseq(pwfup, gene2cat = GOMat, use_genes_without_cat = FALSE) %>%
    as_tibble %>%
    filter(!is.na(ontology))

  termCat <- c('BP', 'MF', 'CC')
  for (j in termCat) {
    write.csv(GOTestWithCatUp %>% filter(ontology == j),
              paste0(i, '_', j, '_up', '.csv') %>% file.path(savepath, .))
  }

  GOTestWithCatDown <- goseq(pwfdown, gene2cat = GOMat, use_genes_without_cat = FALSE) %>%
    as_tibble %>%
    filter(!is.na(ontology))

  termCat <- c('BP', 'MF', 'CC')
  for (j in termCat) {
    write.csv(GOTestWithCatDown %>% filter(ontology == j),
              paste0(i, '_', j, '_down', '.csv') %>% file.path(savepath, .))
    }

  ## KEGG
  pathAnno <- getKEGGPathAnno('ath') %>%
    as_tibble %>%
    mutate(Annotation = Annotation %>% substr(., 1, nchar(.) - 37))

  goseq(pwfup, gene2cat = KEGGMat, use_genes_without_cat = FALSE) %>%
    as_tibble %>%
    inner_join(., pathAnno, by = c('category' = 'pathID')) %>%
    mutate(ontology = 'KEGG') %>%
    write_csv(paste0(i, '_KEGG_up.csv') %>% file.path(savepath, .))

  goseq(pwfdown, gene2cat = KEGGMat, use_genes_without_cat = FALSE) %>%
    as_tibble %>%
    inner_join(., pathAnno, by = c('category' = 'pathID')) %>%
    mutate(ontology = 'KEGG') %>%
    write_csv(paste0(i, '_KEGG_down.csv') %>% file.path(savepath, .))

  ## BioCyc
  pathAnno <- getCycPathway('ARA') %>%
    rename(Annotation = pathAnno) %>%
    mutate(Annotation = Annotation %>% str_replace_all('<.*?>', ''))

  goseq(pwfup, gene2cat = BioCycMat, use_genes_without_cat = FALSE) %>%
    as_tibble %>%
    inner_join(., pathAnno, by = c('category' = 'pathID')) %>%
    mutate(ontology = 'BioCyc') %>%
    write_csv(paste0(i, '_BioCyc_up.csv') %>% file.path(savepath, .))

  goseq(pwfdown, gene2cat = BioCycMat, use_genes_without_cat = FALSE) %>%
    as_tibble %>%
    inner_join(., pathAnno, by = c('category' = 'pathID')) %>%
    mutate(ontology = 'BioCyc') %>%
    write_csv(paste0(i, '_BioCyc_down.csv') %>% file.path(savepath, .))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##############################################################

###################################plot###########################
setwd('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/removeZero/geneset_1stadd/DEGs')

library('ComplexHeatmap')
library('readr')
library('dplyr')
library('magrittr')
library('foreach')
library('RColorBrewer')

cond <- read_csv('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/removeZero/eachGroup_vs_Mock_k_1stadd.csv') %>%
  select(ends_with('padj')) %>%
  colnames %>%
  substring(first = 1, last = nchar(.) - 5)

pathName <- c('BP', 'MF', 'CC', 'KEGG', 'BioCyc')

## up
gsRes <- foreach (i = seq_along(cond), .combine = inner_join) %do% {
  paste0(cond[i], '_', 'BP_up.csv') %>%
    read_csv %>%
    select(2, 7, 6, 3, 5) %>%
    rename(!!paste0(cond[i], '_pvalue') := over_represented_pvalue,
           !!paste0(cond[i], '_in') := numDEInCat,
           size = numInCat)
}

gsResP <- gsRes %>%
  filter(size > 5) %>%
  select(ends_with('pvalue')) %>%
  mutate_all(~ifelse(. > 1, 1, .)) %>%
  mutate_all(~ -log2(.)) %>%
  mutate_all(~ifelse(is.infinite(.), -log2(1e-10), .)) %>%
  mutate_all(~ifelse(. < -log2(0.05), 0, .)) %>% ## no sig --> 0
  setNames(., substring(colnames(.), first = 1, last = nchar(colnames(.)) - 7)) %>%
  {
    sigIdx <- apply(., 1, function(x) {any(x > -log2(0.05))}) %>%
      which ## as least 1 sig
    slice(., sigIdx) %>%
      bind_cols(gsRes %>%
              select(category : size) %>%
              slice(sigIdx))
  }

ht_list <- Heatmap(matrix = gsResP %>% select(-category : -size),
                   name = 'BP',
                   cluster_columns = FALSE,
                   ## row_order = order(scaleC$cl) %>% rev,
                   ## row_split = scaleC$cl,
                   col = colorRampPalette(brewer.pal(n = 5, name = 'Oranges') %>% .[-1:-2] %>% c('white', .))(100),
                   column_title_gp = gpar(fontsize = 7),
                   column_names_gp = gpar(fontsize = 5),
                   use_raster = FALSE)

filePrefix <- 'DEGs_BP_up'

pdf(paste0(filePrefix, '.pdf'))
draw(ht_list)
dev.off()

system(paste0('convert -density 1200 ', paste0(filePrefix, '.pdf'), ' ', paste0(filePrefix, '.jpg')))

write_csv(gsResP, 'DEGs_BP_up')

## down
gsRes <- foreach (i = seq_along(cond), .combine = inner_join) %do% {
  paste0(cond[i], '_', 'BP_down.csv') %>%
    read_csv %>%
    select(2, 7, 6, 3, 5) %>%
    rename(!!paste0(cond[i], '_pvalue') := over_represented_pvalue,
           !!paste0(cond[i], '_in') := numDEInCat,
           size = numInCat)
}

gsResP <- gsRes %>%
  filter(size > 5) %>%
  select(ends_with('pvalue')) %>%
  mutate_all(~ifelse(. > 1, 1, .)) %>%
  mutate_all(~ -log2(.)) %>%
  mutate_all(~ifelse(is.infinite(.), -log2(1e-10), .)) %>%
  mutate_all(~ifelse(. < -log2(0.05), 0, .)) %>%
  setNames(., substring(colnames(.), first = 1, last = nchar(colnames(.)) - 7)) %>%
  {
    sigIdx <- apply(., 1, function(x) {any(x > -log2(0.05))}) %>%
      which ## as least 1 sig
    slice(., sigIdx) %>%
      bind_cols(gsRes %>%
              select(category : size) %>%
              slice(sigIdx))
  }

ht_list <- Heatmap(matrix = gsResP %>% select(-category : -size),
                   name = 'BP',
                   cluster_columns = FALSE,
                   ## row_order = order(scaleC$cl) %>% rev,
                   ## row_split = scaleC$cl,
                   col = colorRampPalette(brewer.pal(n = 6, name = 'PuBu') %>% .[-1:-3] %>% c('white', .))(100),
                   column_title_gp = gpar(fontsize = 7),
                   column_names_gp = gpar(fontsize = 5),
                   use_raster = FALSE)

filePrefix <- 'DEGs_BP_down'

pdf(paste0(filePrefix, '.pdf'))
draw(ht_list)
dev.off()

system(paste0('convert -density 1200 ', paste0(filePrefix, '.pdf'), ' ', paste0(filePrefix, '.jpg')))

write_csv(gsResP, 'DEGs_BP_down')
##################################################################
