##########################cytoscape GO plot###########################
library('tidyverse')
library('magrittr')

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

GOCytoEdge <- function(cpRes) {
  ## INTPUT: `cpRes` is the clusterProfiler table.
  ## OUTPUT: A tibble of edge matrix.
  ## USAGE: Generate the edge table for Cytoscape.

  require('tidyverse')

  ## step1: remove duplicated terms
  cpRes %<>% {
    slice(., which(!duplicated(.$ID)))
  }

  termNum <- nrow(cpRes)

  ## step2: find intersection
  cpResList <- strsplit(cpRes$geneID, split = '/', fixed = TRUE) %>%
    setNames(cpRes$ID)

  interMat <- combWhole_internal(1 : termNum, 1:termNum) %>%
    set_colnames(c('from', 'to')) %>%
    as.data.frame %>%
    as_tibble

  interMat %<>%
    mutate(interNum = apply(., 1, function(x) {
      eachInterNum <- (cpResList[[x[1]]] %in% cpResList[[x[2]]]) %>%
        sum

      return(eachInterNum)
    })) %>%
  filter(interNum > 0) %>%
  mutate(fromAnno = cpRes$ID[from], toAnno = cpRes$ID[to]) %>%
  rename(SOURCE = from, TARGET = to)

  return(interMat)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


setwd('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/removeZero/geneset_1stadd/clusterbc')

cpBP <- read_csv('kmeans10_1stadd_cp_BP.csv')

GOCytoEdge(cpBP) %>% write_csv('tmp1.csv')
#####################################################################

