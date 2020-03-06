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


###############################cluster profiler#####################
library('org.At.tair.db')
library('clusterProfiler')
library('magrittr')
library('tidyverse')

savepath <- '/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/removeZero/geneset/clusterbc'

setwd(savepath)

kmeansRes <- read_csv('../../kmeans10_sig.csv',
                      col_types = cols(Chromosome = col_character()))
kmeansBkg <- read_csv('../../kmeans10.csv',
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
                   pvalueCutoff=0.01,
                   qvalueCutoff=0.01)


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

kall <- lapply(kmeansRes$cl %>% unique, function(x) {

  eachG <- kmeansRes %>% filter(cl == x) %>% .$ID %>% strsplit(split = '.', fixed = TRUE) %>% sapply('[[', 1) %>% unlist %>% unique

  return(eachG)

}) %>%
  set_names(kmeansRes$cl %>% unique %>% paste0('cluster', .))

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
ggsave('kmeans10_cp_BP_dotplot_10.jpg', width = 13)
ggsave('kmeans10_cp_BP_dotplot_10.pdf', width = 13)

kallGOBP %>%
  as.data.frame %>%
  write_csv('kmeans10_cp_BP.csv')

save(kallGOBP, file = 'kmeans10_cp_BP.RData')

kallKEGG <- compareCluster(geneCluster = kall,
                           fun = 'enrichKEGG',
                           organism = 'ath',
                           pvalueCutoff = 0.05)

emapplot(kallGOBP,
         showCategory = 10,
         pie='count',
         pie_scale=1.5,
         layout='nicely')
ggsave('kmeans10_cp_BP_network_10.jpg', width = 18, height = 15)
ggsave('kmeans10_cp_BP_network_10.pdf', width = 18, height = 15)
#######################################################################

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


##############################Plot GO heatmap##########################
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~useful funcs~~~~~~~~~~~~~~~~~~~~~~~~~~~
meanFlg22 <- function(v) {

  require('magrittr')

  res <- v %>%
    split(rep(1 : 4, each = 3)) %>%
    sapply(mean, na.rm = TRUE)

  return(res)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library('tidyverse')
library('DESeq2')
library('ComplexHeatmap')
library('RColorBrewer')
library('circlize')

savepath <- '/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/removeZero/geneset/clusterbc'
setwd(savepath)

load('kmeans10_cp_BP.RData')

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
load('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/removeZero/degres_condi_Mock.RData')

wholeDEG <- read_csv('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/removeZero/eachGroup_vs_Mock_k.csv')
kmeansRes <- read_csv('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/removeZero/kmeans10.csv') %>%
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
interesGO <- list(root_dev = c('GO:0010054', 'GO:0010015', 'GO:0010053', 'GO:0090627'),
                  toxic = c('GO:0009636', 'GO:0098754', 'GO:0010583', 'GO:0009407', 'GO:0009404'),
                  nitrate = c('GO:0015698', 'GO:0071496', 'GO:0006826', 'GO:0046916', 'GO:0015706', 'GO:0010167', 'GO:0000041', 'GO:0031668', 'GO:0098771'),
                  iron = c('GO:0071281', 'GO:0071248', 'GO:0071732', 'GO:0071731', 'GO:0071241'))

## top 10
interesGO <- list(root_dev = c('GO:0010054', 'GO:0010015', 'GO:0010053', 'GO:0090627', 'GO:0048765', 'GO:0048764'),
                  defense = c('GO:0050776', 'GO:0008219', 'GO:0002682', 'GO:0045088'),
                  toxic = c('GO:0009636', 'GO:0098754', 'GO:0010583', 'GO:0009407', 'GO:0009404'),
                  nitrate = c('GO:0015698', 'GO:0006875', 'GO:0009991', 'GO:0030003', 'GO:0046916', 'GO:0098771', 'GO:0000041', 'GO:0015706', 'GO:0019725', 'GO:0010167', 'GO:0055076', 'GO:0006826', 'GO:0031668', 'GO:0010106', 'GO:0055082', 'GO:0071496'),
                  iron = c('GO:1902170', 'GO:1901699', 'GO:0071281', 'GO:0071732', 'GO:0097366', 'GO:0010039', 'GO:0071241', 'GO:0071731', 'GO:0071248', 'GO:0034614'),
                  cell_wall = c('GO:0044036', 'GO:0010410', 'GO:0010383'))


for (i in seq_along(interesGO)) {

  interesGene <- cpBP %>%
    filter(ID %in% interesGO[[i]]) %>%
    .$geneID %>%
    strsplit(split = '/', fixed = TRUE) %>%
    unlist %>%
    unique

  interesMat <- scaleC %>%
    dplyr::filter(GeneID %in% interesGene) %>%
    dplyr::filter(!(cl %in% c(6)))

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
                     column_order = 1 : 4,
                     column_split = rep(c('Mock', 'Mock+flg22', 'Non-sup', 'Sup'), c(1, 1, 1, 1)),
                     show_column_names = FALSE,
                     col = matcol,
                     use_raster = FALSE)

  pdf(paste0(savepath, '/', 'GO_', names(interesGO)[i], '_1stadd_', topGONum, '.pdf'))
  draw(ht_list)
  dev.off()
}
#######################################################################
