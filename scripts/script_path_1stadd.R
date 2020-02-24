#######################GO analysis############################
setwd('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/geneset')

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

kmeansRes <- read_csv('../results/kmeans10_1stadd_sig.csv',
                      col_types = cols(Chromosome = col_character()))
kmeansBkg <- read_csv('../results/kmeans10_1stadd.csv',
                      col_types = cols(Chromosome = col_character()))

savepath <- '/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/removeZero/geneset_1stadd/clusterbc'

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
for (i in kmeansRes$cl %>% unique) {

  prefix <- 'kmeans10_1stadd'

  degVec <- (kmeansRes$cl == i) %>%
    as.integer %>%
    set_names(kmeansRes$ID)

  pwf <- nullp(degVec, bias.data = kmeansRes$Length)

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
    rename(Annotation = pathAnno) %>%
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

  prefix <- 'kmeans10_1stadd'

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
    rename(Annotation = pathAnno) %>%
    mutate(Annotation = Annotation %>% str_replace_all('<.*?>', ''))

  BioCycTestWithCat <- goseq(pwf, gene2cat = BioCycMat, use_genes_without_cat = FALSE) %>%
    as_tibble %>%
    inner_join(., pathAnno, by = c('category' = 'pathID')) %>%
    mutate(ontology = 'BioCyc')

  write.csv(BioCycTestWithCat,
            paste0(prefix, '_cluster', i, '_BioCyc.csv') %>% file.path(savepath, .))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
######################################################################

###################################plot###########################
setwd('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/removeZero/geneset_1stadd/fullbc')

library('ComplexHeatmap')
library('foreach')
library('RColorBrewer')
library('tidyverse')

cond <- read_csv('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/removeZero/kmeans10_1stadd.csv',
                 col_types = cols(Chromosome = col_character())) %>%
  .$cl %>%
  unique %>%
  sort %>%
  paste0('cluster', .)

pathName <- c('BP', 'MF', 'CC', 'KEGG', 'BioCyc')

gsRes <- foreach (i = seq_along(cond), .combine = full_join) %do% {
  paste0('kmeans10_1stadd_', cond[i], '_', 'BP.csv') %>%
    read_csv %>%
    select(2, 7, 6, 3, 5) %>%
    rename(!!paste0(cond[i], '_pvalue') := over_represented_pvalue,
           !!paste0(cond[i], '_in') := numDEInCat,
           size = numInCat)
} %>%
  filter(size > 5) ## filter geneset < 5

gsResP <- gsRes %>%
  select(ends_with('pvalue')) %>%
  mutate_all(~ifelse(is.na(.), 1, .)) %>% ## remove NA
  mutate_all(~ifelse(. > 1, 1, .)) %>% ## remove >1
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
                   col = colorRampPalette(brewer.pal(n = 5, name = 'Oranges') %>% .[-1:-2] %>% c('white', .))(100),
                   column_title_gp = gpar(fontsize = 7),
                   column_names_gp = gpar(fontsize = 5),
                   use_raster = FALSE)

filePrefix <- 'kmeans10_1stadd_BP'

pdf(paste0(filePrefix, '.pdf'))
draw(ht_list)
dev.off()

system(paste0('convert -density 1200 ', paste0(filePrefix, '.pdf'), ' ', paste0(filePrefix, '.jpg')))

write_csv(gsResP, 'kmeans10_1stadd_BP.csv')
####################################################################

###############################cluster profiler#####################
library('org.At.tair.db')
library('clusterProfiler')
library('magrittr')
library('tidyverse')
library('RColorBrewer')

savepath <- '/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/removeZero/geneset_1stadd/clusterbc'

setwd(savepath)

kmeansRes <- read_csv('../../kmeans10_1stadd_sig.csv',
                      col_types = cols(Chromosome = col_character()))
kmeansBkg <- read_csv('../../kmeans10_1stadd.csv',
                      col_types = cols(Chromosome = col_character()))
prefix <- 'kmeans10'

for (i in kmeansRes$cl %>% unique) {

  ## BP
  goBP <- enrichGO(gene = kmeansRes %>% filter(cl == i) %>% .$ID %>% strsplit(split = '.', fixed = TRUE) %>% sapply('[[', 1) %>% unlist %>% unique,
                   OrgDb = 'org.At.tair.db',
                   keyType= 'TAIR',
                   ont = 'BP',
                   universe = keys(org.At.tair.db),
                   pAdjustMethod = 'BH',
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1)


  goBPSim <- clusterProfiler::simplify(goBP,
                                       cutoff = 0.5,
                                       by = 'p.adjust',
                                       select_fun = min)
  ## check and plot
  write.csv(as.data.frame(goBPSim),
            paste0(prefix, '_cluster', i, '_cp_BP.csv') %>% file.path(savepath, .))

  ## KEGG
  kk2 <- enrichKEGG(gene = kmeansRes %>% filter(cl == i) %>% .$ID %>% strsplit(split = '.', fixed = TRUE) %>% sapply('[[', 1) %>% unlist %>% unique,
                    organism = 'ath',
                    pvalueCutoff = 0.05)

  write.csv(as.data.frame(kk2),
            paste0(prefix, '_cluster', i, '_cp_KEGG.csv') %>% file.path(savepath, .))
}

kall <- lapply(kmeansRes$cl %>% unique %>% .[!(. %in% c(9, 10))], function(x) {

  eachG <- kmeansRes %>% filter(cl == x) %>% .$ID %>% strsplit(split = '.', fixed = TRUE) %>% sapply('[[', 1) %>% unlist %>% unique

  return(eachG)

}) %>%
  set_names(kmeansRes$cl %>% unique %>% .[!(. %in% c(9, 10))] %>% paste0('cluster', .))

kallGOBP <- compareCluster(geneCluster = kall,
                           fun = 'enrichGO',
                           OrgDb = 'org.At.tair.db',
                           keyType= 'TAIR',
                           ont = 'BP',
                           universe = keys(org.At.tair.db),
                           pAdjustMethod = 'BH',
                           pvalueCutoff=0.01,
                           qvalueCutoff=0.1)

kallGOBPSim <- clusterProfiler::simplify(kallGOBP,
                                         cutoff = 0.9,
                                         by = 'p.adjust',
                                         select_fun = min)
dotplot(kallGOBPSim, showCategory = 20)

dotplot(kallGOBP)
ggsave('kmeans10_1stadd_cp_BP_dotplot.jpg', width = 13)
ggsave('kmeans10_1stadd_cp_BP_dotplot.pdf', width = 13)

kallGOBP %>%
  as.data.frame %>%
  write_csv('kmeans10_1stadd_cp_BP.csv')

save(kallGOBP, file = 'kmeans10_1stadd_cp_BP.RData')

kallGOBPPlot <- kallGOBP
kallGOBPPlot@compareClusterResult %<>% {
  clusterColor <- c(brewer.pal(n = 8, name = 'Set1'))[.$Cluster]
  cbind(., clusterColor)
}

emapplot(kallGOBP,
         showCategory = 5,
         pie='count',
         pie_scale=1.5,
         layout='nicely')
ggsave('kmeans10_1stadd_cp_BP_network.jpg', width = 18, height = 15)
ggsave('kmeans10_1stadd_cp_BP_network.pdf', width = 18, height = 15)

kallKEGG <- compareCluster(geneCluster = kall,
                           fun = 'enrichKEGG',
                           organism = 'ath',
                           pvalueCutoff = 0.05)
dotplot(kallKEGG)
#######################################################################

###############################metacape################################
library('foreach')
library('tidyverse')
library('magrittr')

setwd('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/removeZero/geneset_1stadd/clusterbc/')

kcluster <- read_csv('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/removeZero/kmeans10_1stadd_sig.csv')

clusterSeq <- kcluster$cl %>% unique %>% seq_along
kall <- foreach (i = clusterSeq) %do% {

  eachTrans <- kcluster %>%
    filter(cl == i) %>%
    .$ID

  return(eachTrans)
}

maxLen <- sapply(kall, length) %>% max

kall %<>% lapply(function(x) {
  x %<>% c(., rep('', maxLen - length(x)))
  return(x)
}) %>%
  do.call(cbind, .) %>%
  set_colnames(paste0('cluster', clusterSeq))

write_csv(kall %>% as_tibble %>% select(-cluster9:-cluster10), 'kmeans10_1stadd_1to8_sig_metascape.csv')
write_csv(kcluster %>% select(ID), 'kall_bkg.csv')
#######################################################################
