##########################cytoscape GO plot###########################
library('tidyverse')
library('reshape2')
library('magrittr')
library('foreach')

##~~~~~~~~~~~~~~~~~~~~~useful funcs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
comb2_internal <- function(x) {

  l <- length(x)

  if (l < 2) {
    m <- matrix(0, nrow = 0, ncol = 2)
  } else {
    fV <- rep(x[1:(l - 1)], ((l - 1):1))
    tVIdx <- unlist(mapply(seq, 2:l, l))
    tV <- x[tVIdx]

    m <- cbind(fV, tV)
  }

  return(m)
}

combWhole_internal <- function(x, y, self = FALSE, bidirect = FALSE) {

  ## check x
  ## first unique x
  x <- unique(x)
  x <- x[x %in% y]

  l <- length(x)

  if (l >= length(y)) {
    m <- comb2_internal(x)
  } else {
    yLeft <- y[!(y %in% x)]
    m <- cbind(rep(x, each = length(yLeft)),
               rep(yLeft, l))
    m <- rbind(m,
               comb2_internal(x))
  }

  if (bidirect) {
    m <- rbind(m, m[, 2:1])
  } else {}

  if (self) {
    m <- rbind(m, combSelf_internal(x))
  } else {}

  colnames(m) <- NULL

  return(m)
}

JacSim <- function (x, y) {
  x <- unlist(x)
  y <- unlist(y)
  length(intersect(x, y))/length(unique(c(x, y)))
}

Shrinkage <- function(x, minNew, maxNew) {
  ## INPUT: `x` is numeric vector. `minNew` and `maxNew` is the new min and new max.
  ## OUTPUT: A shrink vector.

  unit <- (maxNew - minNew) / (max(x) - min(x))

  res <- minNew + (x - min(x)) * unit

  return(res)
}

GOCytoEdge <- function(cpRes, JacSimThres = 0.2) {
  ## INTPUT: `cpRes` is the clusterProfiler table. `JacSimThres` is the threshold of Jaccard simialrity.
  ## OUTPUT: A tibble of edge matrix.
  ## USAGE: Generate the edge table for Cytoscape.

  require('tidyverse')

  ## step1: remove duplicated terms
  cpRes %<>%
    select(ID, geneID, Description) %>%
    group_by(ID) %>%
    summarise(geneID = paste(geneID, collapse = '/'),
              Description = sample(Description, 1)) %>%
    ungroup

  termNum <- nrow(cpRes)

  ## step2: find intersection
  cpResList <- strsplit(cpRes$geneID, split = '/', fixed = TRUE) %>%
    sapply(unique) %>%
    setNames(cpRes$ID)

  interMat <- combWhole_internal(1 : termNum, 1:termNum) %>%
    set_colnames(c('from', 'to')) %>%
    as.data.frame %>%
    as_tibble

  interMat %<>%
    mutate(jacSim = apply(., 1, function(x) {
      eachJacSim <- JacSim(cpResList[[x[1]]], cpResList[[x[2]]])
      return(eachJacSim)
    })) %>%
    filter(jacSim >= JacSimThres) %>% ## filter by jaccard similarity
    mutate(edgeWidth = Shrinkage(jacSim, 2, 5)) %>% ## edge width
    mutate(SOURCE = cpRes$ID[from], TARGET = cpRes$ID[to]) %>%
    mutate(fromDesc = cpRes$Description[from], toDesc = cpRes$Description[to]) %>%
    select(-from, -to)

  return(interMat)
}

GOCytoNode <- function(cpRes) {
  ## INTPUT: `cpRes` is the clusterProfiler table.
  ## OUTPUT: A tibble of node matrix.
  ## USAGE: Generate the node table for Cytoscape.

  require('reshape2')

  ## unique GO terms
  cpUniq <- cpRes %>%
    select(ID, geneID, geneName, Description) %>%
    group_by(ID) %>%
    summarise(geneID = paste(geneID, collapse = '/'),
              Description = sample(Description, 1),
              geneID = sample(geneID, 1),
              geneName = sample(geneName, 1)) %>%
    ungroup

  nodeMat <- cpRes %>%
    mutate(Cluster = str_extract(Cluster, 'cluster\\d')) %>%
    dcast(ID ~ Cluster, value.var = 'Count') %>%
    mutate_all(~ifelse(is.na(.), 0, .)) %>%
    mutate(nodeSize = select(., -ID) %>% rowSums %>% Shrinkage(30, 50)) %>%
    inner_join(cpUniq) %>%
    rename(SOURCE=ID)

  return(nodeMat)

}

CytoGeneEdge_ <- function(GOGene, geneAnno) {

  ## INPUT: `GOgene` is the 'GO-gene' matrix. `geneAnno` is the gene annotation.
  ## OUTPUT: A `tibble` of interaction matrix.

  res <- tibble(SOURCE = strsplit(GOGene$geneName, split = '/', fixed = TRUE) %>% unlist,
                TARGET = GOGene$ID,
                toDesc = GOGene$Description) %>%
    mutate(jacSim = 1,
           edgeWidth = 1) %>%
    mutate(fromDesc = anno$Description[match(SOURCE, anno$Gene)]) %>%
    select(jacSim, edgeWidth, SOURCE, TARGET, fromDesc, toDesc)

  return(res)
}


GOCytoGeneEdge <- function(cpRes, ...) {
  ## INPUT: `cpRes` is the clusterProfiler table.
  ## OUTPUT: A tibble of gene-GO matrix.
  ## USAGE: Generate the node table for Cytoscape.

  require('foreach')

  res <- foreach(i = seq_len(nrow(cpRes)), .combine = bind_rows) %do% {
    CytoGeneEdge_(cpRes[i, ], ...)
  }

  return(res)
}

CytoGeneNode_ <- function(GOGene, anno) {
  ## INPUT: `GOgene` is the 'GO-gene' matrix. `geneAnno` is the gene annotation.
  ## OUTPUT: A `tibble` of node matrix.

  GOGene %<>%
    mutate(Cluster = str_extract(Cluster, 'cluster\\d'))

  res <- tibble(SOURCE = strsplit(GOGene$geneName, split = '/', fixed = TRUE) %>% unlist,
                Cluster = GOGene$Cluster) %>%
    mutate(nodeSize = 10,
      geneID = strsplit(GOGene$geneID, split = '/', fixed = TRUE) %>% unlist,
      Description = anno$Description[match(SOURCE, anno$Gene)],
      geneName = SOURCE)

  return(res)
}


GOCytoGeneNode <- function(cpRes, ...) {
  ## INPUT: `cpRes` is the clusterProfiler table.
  ## OUTPUT: A tibble of gene-GO edge.
  ## USAGE: Generate the node table for Cytoscape.

  require('foreach')
  require('reshape2')

  nodeMat <- foreach(i = seq_len(nrow(cpRes)), .combine = bind_rows) %do% {
    CytoGeneNode_(cpRes[i, ], ...)
  } %>%
    slice(which(!duplicated(.)))

  clusterMat <- nodeMat %>%
    mutate(Val = 1) %>%
    dcast(SOURCE ~ Cluster, value.var = 'Val', fill = 0) %>%
    mutate_at(.vars = vars(contains('cluster')),
              .funs = ~ifelse(. > 0, 1, .))

  res <- inner_join(clusterMat, nodeMat %>% select(-Cluster)) %>%
    as_tibble

  return(res)
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/removeZero/geneset_1stadd/clusterbc')

load('kmeans10_1stadd_cp_BP.RData')

anno <- read_csv('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/Ensembl_ath_Anno.csv',
                 col_types = cols(Chromosome = col_character())) %>%
  mutate_all(list(~replace(., is.na(.), ''))) %>%
  mutate(GeneID = strsplit(ID, split = '.', fixed = TRUE) %>%
           sapply('[[', 1) %>%
           unlist) %>%
  mutate(Gene = if_else(nchar(Gene) == 0, GeneID, Gene)) %>%
  select(GeneID, Gene, Description) %>%
  slice(which(!duplicated(.)))


cpBP <- clusterProfiler:::fortify.compareClusterResult(kallGOBP,
                                                       showCategory = 10) %>%
  as_tibble %>%
  mutate(geneName = sapply(geneID, function(x) {
    strsplit(x, split = '/', fixed = TRUE) %>%
      unlist %>%
      tibble(GeneID = .) %>%
      inner_join(anno) %>%
      .$Gene %>%
      paste(collapse = '/')
  }))

## GO network
GOCytoEdge(cpBP) %>% write_csv('tmp1.csv')
GOCytoNode(cpBP) %>% write_csv('tmp2.csv')
cpBP %>% write_csv('tmp3.csv')

## GO-gene network
bind_rows(GOCytoEdge(cpBP),
          GOCytoGeneEdge(cpBP, anno)) %>% write_csv('tmp1.csv')

bind_rows(GOCytoNode(cpBP),
          GOCytoGeneNode(cpBP, anno)) %>% write_csv('tmp2.csv')

##~~~~~~~~~~~~~~~~interesting genes cluster~~~~~~~~~~~~~~~~~~~~~~~~~~
## AT2G19190 -- SIRK/FRK1 -- 3
## AT1G14550 -- PER5 -- 3
## AT1G18570 -- MYB51 -- 5
## AT2G30750 -- CYP71A12 -- 5
## AT1G19250 -- FMO1 -- 5
## AT1G73805 -- SARD1 -- 1
## AT4G28460 -- PIP1 -- 5
## AT4G37290 -- PIP2 -- 5
## AT3G48090 -- EDS1 -- 1
## AT3G07040 -- RPM1 -- 1
## AT3G50950 -- ZAR1/RPP13L4 -- 1
## AT4G11170 -- -- 5
## AT5G41540 -- -- 1
interesGene <- c('AT2G19190',
                 'AT1G14550',
                 'AT1G18570',
                 'AT2G30750',
                 'AT1G73805',
                 'AT4G28460',
                 'AT4G37290',
                 'AT3G48090',
                 'AT3G07040',
                 'AT3G50950',
                 'AT4G11170',
                 'AT5G41540')

str_detect(GOCytoNode(cpBP)$geneID, interesGene %>% paste(collapse = '|')) %>%
  GOCytoNode(cpBP)$Description[.]
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#####################################################################

