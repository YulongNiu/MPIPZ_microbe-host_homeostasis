
## originally by Yulong Niu
## yulong.niu@hotmail.com

#######################GO analysis############################
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

kmeansRes <- read_csv('../results/inter_col0mu_35down_1stadd.csv',
                      col_types = cols(Chromosome = col_character())) %>%
  select(ID, Length, cl)

athAnno <- read_csv('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/Ensembl_ath_Anno.csv',
                    col_types = cols(Chromosome = col_character())) %>%
  mutate(Gene = Gene %>% {if_else(is.na(.), '', .)}) %>%
  mutate(Description = Description %>% {if_else(is.na(.), '', .)})

##~~~~~~~~~~~~~~~select genesets~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
athGO %<>%
  lapply(function(x) {x[x %in% athAnno$ID]}) %>%
  {
    l <- sapply(., length) > 0
    .[l]
  }
GOMat <- foreach(i = 1:length(athGO), .combine = rbind) %dopar% {
  eachMat <- cbind(athGO[[i]], names(athGO)[i])
  return(eachMat)
} %>% as.data.frame

athKEGG %<>%
  lapply(function(x) {x[x %in% athAnno$ID]}) %>%
  {
    l <- sapply(., length) > 0
    .[l]
  }
KEGGMat <- foreach(i = 1:length(athKEGG), .combine = rbind) %dopar% {
  eachMat <- cbind(athKEGG[[i]], names(athKEGG)[i])
  return(eachMat)
} %>% as.data.frame


athBioCyc %<>%
  lapply(function(x) {x[x %in% athAnno$ID]}) %>%
  {
    l <- sapply(., length) > 0
    .[l]
  }
BioCycMat <- foreach(i = 1:length(athBioCyc), .combine = rbind) %dopar% {
  eachMat <- cbind(athBioCyc[[i]], names(athBioCyc)[i])
  return(eachMat)
} %>% as.data.frame
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~choose sig~~~~~~~~~~~~~~~~~~~~~
sig <- (athAnno$ID %in% kmeansRes$ID) %>%
  as.numeric %>%
  set_names(athAnno$ID)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~GO~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pwf <- nullp(sig, bias.data = athAnno$Length)

GOTestWithCat <- goseq(pwf, gene2cat = GOMat, use_genes_without_cat = FALSE) %>%
  as_tibble %>%
  filter(!is.na(ontology)) %>%
  rename(Annotation = term)

termCat <- c('BP', 'MF', 'CC')
for (k in termCat) {
  write.csv(GOTestWithCat %>% filter(ontology == k),
            paste0('kmeans10_', k, '_35up_all.csv') %>% file.path('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/pathway_35down_all', .))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~KEGG~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pathAnno <- getKEGGPathAnno('ath') %>%
  as_tibble %>%
  mutate(Annotation = Annotation %>% substr(., 1, nchar(.) - 37))

pwf <- nullp(sig, bias.data = kmeansRes$Length)

KEGGTestWithCat <- goseq(pwf, gene2cat = KEGGMat, use_genes_without_cat = FALSE) %>%
  as_tibble %>%
  inner_join(., pathAnno, by = c('category' = 'pathID')) %>%
  mutate(ontology = 'KEGG')

write.csv(KEGGTestWithCat,
          paste0('kmeans10_KEGG_35up_all.csv') %>% file.path('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/pathway_35down_all', .))
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~BioCyc~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pathAnno <- getCycPathway('ARA') %>%
  rename(Annotation = pathAnno) %>%
  mutate(Annotation = Annotation %>% str_replace_all('<.*?>', ''))

pwf <- nullp(sig, bias.data = kmeansRes$Length)

BioCycTestWithCat <- goseq(pwf, gene2cat = BioCycMat, use_genes_without_cat = FALSE) %>%
  as_tibble %>%
  inner_join(., pathAnno, by = c('category' = 'pathID')) %>%
  mutate(ontology = 'BioCyc')

write.csv(BioCycTestWithCat,
          paste0('kmeans10_BioCyc_35up_all.csv') %>% file.path('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/pathway_35down_all', .))
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##############################################################


###################################plot###########################
setwd('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/pathway_35up_1stadd')

library('ggplot2')
library('readr')
library('dplyr')
library('magrittr')
library('foreach')

cln <- 1:10
cln <- 10
cln <- 1

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
               numDEInCat >= 1)
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
