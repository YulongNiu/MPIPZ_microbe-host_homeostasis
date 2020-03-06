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

dotplot(kallGOBP, showCategory = 10)
ggsave('kmeans10_1stadd_cp_BP_dotplot_10.jpg', width = 13)
ggsave('kmeans10_1stadd_cp_BP_dotplot_10.pdf', width = 13)

kallGOBP %>%
  as.data.frame %>%
  write_csv('kmeans10_1stadd_cp_BP.csv')

save(kallGOBP, file = 'kmeans10_1stadd_cp_BP.RData')

emapplot(kallGOBP,
         showCategory = 5,
         pie='count',
         pie_scale=1.5,
         layout='nicely')
ggsave('kmeans10_1stadd_cp_BP_network.jpg', width = 18, height = 16)
ggsave('kmeans10_1stadd_cp_BP_network.pdf', width = 18, height = 16)

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

##############################Plot GO heatmap##########################
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~useful funcs~~~~~~~~~~~~~~~~~~~~~~~~~~~
meanFlg22 <- function(v) {

  require('magrittr')

  res <- v %>%
    split(rep(1 : 10, each = 4)) %>%
    sapply(mean, na.rm = TRUE)

  return(res)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library('tidyverse')
library('DESeq2')
library('ComplexHeatmap')
library('RColorBrewer')
library('circlize')

savepath <- '/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/removeZero/geneset_1stadd/clusterbc'
setwd(savepath)

load('kmeans10_1stadd_cp_BP.RData')

topGONum <- 10

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Sig terms~~~~~~~~~~~~~~~~~~~~~~~~~~~~
anno <- read_csv('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/Ensembl_ath_Anno.csv',
                 col_types = cols(Chromosome = col_character())) %>%
  mutate_all(list(~replace(., is.na(.), ''))) %>%
  mutate(GeneID = strsplit(ID, split = '.', fixed = TRUE) %>%
           sapply('[[', 1) %>%
           unlist) %>%
  mutate(Gene = if_else(nchar(Gene) == 0, GeneID, Gene)) %>%
  dplyr::select(GeneID, Gene, Description) %>%
  dplyr::slice(which(!duplicated(.)))

cpBP <- clusterProfiler:::fortify.compareClusterResult(kallGOBP,
                                                       showCategory = topGONum) %>%
  as_tibble %>%
  mutate(geneName = sapply(geneID, function(x) {
    strsplit(x, split = '/', fixed = TRUE) %>%
      unlist %>%
      tibble(GeneID = .) %>%
      inner_join(anno) %>%
      .$Gene %>%
      paste(collapse = '/') %>%
      str_replace('C/VIF2', 'C-VIF2') ## replace genes with '/'
  }))
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~heatmap matrix~~~~~~~~~~~~~~~~~~~~~~
load('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/removeZero/degres_condi_Mock_1stadd.RData')

wholeDEG <- read_csv('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/removeZero/eachGroup_vs_Mock_k_1stadd.csv')
kmeansRes <- read_csv('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/removeZero/kmeans10_1stadd.csv') %>%
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

rawC <- rldData %>%
  as.data.frame %>%
  rownames_to_column('ID') %>%
  as_tibble %>%
  dplyr::select(matches('Mock_\\d|HKSynCom33_\\d|HKSynCom35_\\d'), matches('Mock_Flg22_\\d|HKSynCom33_Flg22_\\d|HKSynCom35_Flg22_\\d'), everything()) %>%
  inner_join(heatsig %>% select(ID, cl))
## inner_join(kmeansRes) ## all transcripts

## scale counts
scaleC <- rawC %>%
  select(contains('_')) %>%
  t %>%
  scale %>%
  t %>%
  as_tibble %>%
  bind_cols(rawC %>% select(ID, cl)) %>%
  mutate(GeneID = ID %>% strsplit(split = '.', fixed = TRUE) %>% sapply('[', 1))
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## top 5
interesGO <- list(root_dev = c('GO:0010053', 'GO:0010054', 'GO:0010015', 'GO:0048764', 'GO:0090627'),
                  defense = c('GO:0045730', 'GO:0002679', 'GO:0010200', 'GO:0010243', 'GO:0009753'),
                  hypoxia = c('GO:0036294', 'GO:0071453', 'GO:0071456'),
                  toxic = c('GO:0009636', 'GO:0098754', 'GO:0010583', 'GO:0009407', 'GO:0009404'),
                  nitrate = c('GO:0015706', 'GO:0010167', 'GO:0015698', 'GO:0000041'),
                  cell_wall = c('GO:0045491', 'GO:0045492', 'GO:0042546', 'GO:0044038', 'GO:0010413'))

## top 10
interesGO <- list(root_dev = c('GO:0090627', 'GO:0080147', 'GO:0048767', 'GO:0048765', 'GO:0048764', 'GO:0048588', 'GO:0048469', 'GO:0010054', 'GO:0010053', 'GO:0010015'),
                  defense = c('GO:0050832', 'GO:0046189', 'GO:0045730', 'GO:0035690', 'GO:0010243', 'GO:0010200', 'GO:0009753', 'GO:0009723', 'GO:0009697', 'GO:0009611', 'GO:0002679'),
                  hypoxia = c('GO:0036294', 'GO:0071453', 'GO:0071456'),
                  toxic = c('GO:0009636', 'GO:0098754', 'GO:0010583', 'GO:0009407', 'GO:0009404'),
                  nitrate = c('GO:0048878', 'GO:0042594', 'GO:0031669', 'GO:0030003', 'GO:0019725', 'GO:0015706', 'GO:0015698', 'GO:0010167', 'GO:0009267', 'GO:0006875', 'GO:0006826', 'GO:0000041'),
                  cell_wall = c('GO:0070592', 'GO:0070589', 'GO:0045492', 'GO:0045491', 'GO:0044038', 'GO:0044036', 'GO:0042546', 'GO:0010413', 'GO:0010410', 'GO:0010383'))


for (i in seq_along(interesGO)) {

  interesGene <- cpBP %>%
    filter(ID %in% interesGO[[i]]) %>%
    .$geneID %>%
    strsplit(split = '/', fixed = TRUE) %>%
    unlist %>%
    unique

  interesMat <- scaleC %>%
    dplyr::filter(GeneID %in% interesGene) %>%
    dplyr::filter(!(cl %in% c(9:10)))

  matcol <- colorRamp2(seq(min(scaleC %>% select(contains('_'))), max(scaleC %>% select(contains('_'))), length = 100), colorRampPalette(rev(brewer.pal(n = 10, name = 'Spectral'))[c(-3, -4, -6, -7)])(100))

  dim(interesMat) %>% print

  ht_list <- Heatmap(matrix = interesMat %>%
                       select(contains('_')) %>%
                       apply(1, meanFlg22) %>%
                       t,
                     name = 'Scaled Counts',
                     row_order = order(interesMat$cl) %>% rev,
                     row_split = interesMat$cl,
                     row_gap = unit(2, "mm"),
                     column_order = 1 : 10,
                     column_split = rep(c('Mock/HKSynCom', 'Non-sup', 'Sup'), c(6, 2, 2)),
                     show_column_names = FALSE,
                     col = matcol,
                     use_raster = FALSE)

  pdf(paste0(savepath, '/', 'GO_', names(interesGO)[i], '_1stadd_', topGONum, '.pdf'))
  draw(ht_list)
  dev.off()
}
#######################################################################
