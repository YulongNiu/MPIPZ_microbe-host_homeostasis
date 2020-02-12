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

kmeansRes <- read_csv('../results/removeZero/kmeans10_sig.csv',
                      col_types = cols(Chromosome = col_character()))
kmeansBkg <- read_csv('../results/removeZero/kmeans10.csv',
                      col_types = cols(Chromosome = col_character()))

savepath <- '/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/removeZero/geneset/fullbc/'

##~~~~~~~~~~~~~~~select genesets~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
athGO %<>%
  lapply(function(x) {x[x %in% kmeansBkg$ID]}) %>%
  {
    l <- sapply(., length) > 0
    .[l]
  }
GOMat <- foreach(i = 1:length(athGO), .combine = rbind) %dopar% {
  eachMat <- cbind(athGO[[i]], names(athGO)[i])
  return(eachMat)
} %>% as.data.frame

athKEGG %<>%
  lapply(function(x) {x[x %in% kmeansBkg$ID]}) %>%
  {
    l <- sapply(., length) > 0
    .[l]
  }
KEGGMat <- foreach(i = 1:length(athKEGG), .combine = rbind) %dopar% {
  eachMat <- cbind(athKEGG[[i]], names(athKEGG)[i])
  return(eachMat)
} %>% as.data.frame

athBioCyc %<>%
  lapply(function(x) {x[x %in% kmeansBkg$ID]}) %>%
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
for (i in kmeansBkg$cl %>% unique) {

  prefix <- 'kmeans10'

  degVec <- (kmeansBkg$cl == i) %>%
    as.integer %>%
    set_names(kmeansBkg$ID)

  pwf <- nullp(degVec, bias.data = kmeansBkg$Length)

  ## GO
  GOTestWithCat <- goseq(pwf, gene2cat = GOMat, use_genes_without_cat = FALSE) %>%
    as_tibble %>%
    filter(!is.na(ontology))

  termCat <- c('BP', 'MF', 'CC')
  for (j in termCat) {
    write.csv(GOTestWithCat %>% filter(ontology == j),
              paste0(prefix, '_cluster', i, '_', j, '.csv') %>% file.path(savepath, .))
  }

  ## KEGG
  pathAnno <- getKEGGPathAnno('ath') %>%
    as_tibble %>%
    mutate(Annotation = Annotation %>% substr(., 1, nchar(.) - 37))

  KEGGTestWithCat <- goseq(pwf, gene2cat = KEGGMat, use_genes_without_cat = FALSE) %>%
    as_tibble %>%
    inner_join(., pathAnno, by = c('category' = 'pathID')) %>%
    mutate(ontology = 'KEGG')

  write.csv(KEGGTestWithCat,
            paste0(prefix, '_cluster', i, '_KEGG.csv') %>% file.path(savepath, .))

  ## BioCyc
  pathAnno <- getCycPathway('ARA') %>%
    dplyr::rename(Annotation = pathAnno) %>%
    mutate(Annotation = Annotation %>% str_replace_all('<.*?>', ''))

  BioCycTestWithCat <- goseq(pwf, gene2cat = BioCycMat, use_genes_without_cat = FALSE) %>%
    as_tibble %>%
    inner_join(., pathAnno, by = c('category' = 'pathID')) %>%
    mutate(ontology = 'BioCyc')

  write.csv(BioCycTestWithCat,
            paste0(prefix, '_cluster', i, '_BioCyc.csv') %>% file.path(savepath, .))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~whole cluster gene-set with background~~~~~~~~~~
for (i in kmeansRes$cl %>% unique) {

  prefix <- 'kmeans10'

  eachRes <- kmeansRes %>%
    filter(cl == i) %>%
    {.$ID}
  eachBkg <- kmeansBkg %>%
    filter(cl == i) %>%
    {.$ID}
  eachLength <- kmeansBkg %>%
    filter(cl == i) %>%
    {.$Length}

  degVec <- rep(0, length(eachBkg)) %>%
    set_names(eachBkg)
  degVec[match(eachRes, eachBkg)] <- 1

  pwf <- nullp(degVec, bias.data = eachLength)

  ## GO
  GOTestWithCat <- goseq(pwf, gene2cat = GOMat, use_genes_without_cat = FALSE) %>%
    as_tibble %>%
    filter(!is.na(ontology))

  termCat <- c('BP', 'MF', 'CC')
  for (j in termCat) {
    write.csv(GOTestWithCat %>% filter(ontology == j),
              paste0(prefix, '_cluster', i, '_', j, '.csv') %>% file.path(savepath, .))
  }

  ## KEGG
  pathAnno <- getKEGGPathAnno('ath') %>%
    as_tibble %>%
    mutate(Annotation = Annotation %>% substr(., 1, nchar(.) - 37))

  KEGGTestWithCat <- goseq(pwf, gene2cat = KEGGMat, use_genes_without_cat = FALSE) %>%
    as_tibble %>%
    inner_join(., pathAnno, by = c('category' = 'pathID')) %>%
    mutate(ontology = 'KEGG')

  write.csv(KEGGTestWithCat,
            paste0(prefix, '_cluster', i, '_KEGG.csv') %>% file.path(savepath, .))

  ## BioCyc
  pathAnno <- getCycPathway('ARA') %>%
    dplyr::rename(Annotation = pathAnno) %>%
    mutate(Annotation = Annotation %>% str_replace_all('<.*?>', ''))

  BioCycTestWithCat <- goseq(pwf, gene2cat = BioCycMat, use_genes_without_cat = FALSE) %>%
    as_tibble %>%
    inner_join(., pathAnno, by = c('category' = 'pathID')) %>%
    mutate(ontology = 'BioCyc')

  write.csv(BioCycTestWithCat,
            paste0(prefix, '_cluster', i, '_BioCyc.csv') %>% file.path(savepath, .))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~choose sig~~~~~~~~~~~~~~~~~~~~~
## padj < 0.05 & |log2FC| > log2(1.5)
fcsig <- kmeansRes %>%
  select(ends_with('FoldChange')) %>%
  transmute_all(list(~ case_when(. > 1 ~ 1,
                                 . < -1 ~ -1,
                                 TRUE ~ 0)))
padjsig <- kmeansRes %>%
  select(ends_with('padj')) %>%
  abs %>%
  `<`(0.05) %>%
  as_tibble %>%
  transmute_all(list(~ if_else(is.na(.), FALSE, .)))

sig <- (padjsig * fcsig) %>%
  as_tibble %>%
  mutate(ID = kmeansRes$ID, cl = kmeansRes$cl) %>%
  select(ID, everything()) %>%
  rename(Flg22_vs_Mock = Flg22_vs_Mock_padj,
         Flg22_SynCom33_vs_Mock = Flg22_SynCom33_vs_Mock_padj,
         Flg22_SynCom35_vs_Mock = Flg22_SynCom35_vs_Mock_padj)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~GO~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
vsGroup <- c('Flg22_vs_Mock', 'Flg22_SynCom33_vs_Mock', 'Flg22_SynCom35_vs_Mock')
cln <- 1:10
cln <- 4

for (i in vsGroup) {
  for (j in cln) {
    degVec <- sig %>%
      transmute((!!as.name(i)) != 0 &
                cl == j) %>%
      unlist %>%
      as.integer
    names(degVec) <- sig$ID

    pwf <- nullp(degVec, bias.data = kmeansRes$Length)

    GOTestWithCat <- goseq(pwf, gene2cat = GOMat, use_genes_without_cat = FALSE) %>%
      as_tibble %>%
      filter(!is.na(ontology)) %>%
      rename(Annotation = term)

    termCat <- c('BP', 'MF', 'CC')
    for (k in termCat) {
      write.csv(GOTestWithCat %>% filter(ontology == k),
                paste0('kmeans10_', i, '_cluster', j, '_', k, '.csv') %>% file.path('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/pathway_35up', .))
    }
  }
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~KEGG~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pathAnno <- getKEGGPathAnno('ath') %>%
  as_tibble %>%
  mutate(Annotation = Annotation %>% substr(., 1, nchar(.) - 37))

for (i in vsGroup) {
  for (j in cln) {
    degVec <- sig %>%
      transmute((!!as.name(i)) != 0 &
                cl == j) %>%
      unlist %>%
      as.integer
    names(degVec) <- sig$ID

    pwf <- nullp(degVec, bias.data = kmeansRes$Length)

    KEGGTestWithCat <- goseq(pwf, gene2cat = KEGGMat, use_genes_without_cat = FALSE) %>%
      as_tibble %>%
      inner_join(., pathAnno, by = c('category' = 'pathID')) %>%
      mutate(ontology = 'KEGG')

    write.csv(KEGGTestWithCat,
              paste0('kmeans10_', i, '_cluster', j, '_KEGG', '.csv') %>% file.path('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/pathway_35up', .))
  }
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~BioCyc~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pathAnno <- getCycPathway('ARA') %>%
  rename(Annotation = pathAnno) %>%
  mutate(Annotation = Annotation %>% str_replace_all('<.*?>', ''))

for (i in vsGroup) {
  for (j in cln) {
    degVec <- sig %>%
      transmute((!!as.name(i)) != 0 &
                cl == j) %>%
      unlist %>%
      as.integer
    names(degVec) <- sig$ID

    pwf <- nullp(degVec, bias.data = kmeansRes$Length)

    BioCycTestWithCat <- goseq(pwf, gene2cat = BioCycMat, use_genes_without_cat = FALSE) %>%
      as_tibble %>%
      inner_join(., pathAnno, by = c('category' = 'pathID')) %>%
      mutate(ontology = 'BioCyc')

    write.csv(BioCycTestWithCat,
              paste0('kmeans10_', i, '_cluster', j, '_BioCyc', '.csv') %>% file.path('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/pathway_35up', .))
  }
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##############################################################


###################################plot###########################
setwd('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/pathway_35up')

library('ggplot2')
library('readr')
library('dplyr')
library('magrittr')
library('foreach')

cln <- 1:10
cln <- 4
geneset <- c('BP', 'MF', 'CC', 'KEGG', 'BioCyc')

for (i in cln) {
  for (j in geneset) {

    vsGroup <- c('Flg22_vs_Mock', 'Flg22_SynCom33_vs_Mock', 'Flg22_SynCom35_vs_Mock')

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
