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

library('ggplot2')
library('readr')
library('dplyr')
library('magrittr')
library('foreach')

## cln <- 1:10
## cln <- 10
## cln <- 1
cln <- 2

geneset <- c('BP', 'MF', 'CC', 'KEGG', 'BioCyc')

for (i in cln) {
  for (j in geneset) {

    vsGroup <- c('Mock_Flg22_vs_Mock', 'SynCom33_Flg22_vs_Mock', 'SynCom35_Flg22_vs_Mock')

    pathPlot <- foreach(k = seq_along(vsGroup), .combine = bind_rows) %do% {
      vsGroup[k] %>%
        {paste0('kmeans10_', ., '_cluster', i, '_', j, '.csv')} %>%
        read_csv %>%
        select(Annotation, over_represented_pvalue, numDEInCat, numInCat) %>%
        rename(pvalue = over_represented_pvalue) %>%
        mutate(group = vsGroup[k], ratio = numDEInCat / numInCat) %>%
        filter(pvalue < 0.05 &
               numDEInCat >= 2)
    }

    colorPal <- colorRampPalette(rev(c('red', 'yellow', 'cyan', 'blue')), bias=1)(10)

    ggplot(pathPlot, aes(x = group, y = Annotation)) +
      geom_point(aes(size = ratio, colour = -log10(pvalue))) +
      scale_colour_gradientn(name = '-log10(P-value)', limits=c(0, max(-log10(pathPlot$pvalue))), colours = colorPal) +
      ## scale_x_discrete(labels = c('Flg22', 'Flg22+SynCom33', 'Flg22+SynCom35')) +
      ylab(j) +
      xlab('') +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
    ggsave(paste0('kmeans10_cluster', i, '_', j, '.pdf'))
    ggsave(paste0('kmeans10_cluster', i, '_', j, '.jpg'))
  }
}
##################################################################
