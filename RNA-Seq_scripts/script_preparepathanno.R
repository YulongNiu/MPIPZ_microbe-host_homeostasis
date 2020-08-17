
## originally by Yulong Niu
## yulong.niu@hotmail.com

#########################Prepare GO from Consortium###################
library('readr')
library('dplyr')
library('stringr')

setwd('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/')

athGORaw <- read_delim('ath_GO_raw.txt',
                       delim = '\t',
                       col_names = c('TAIRLoc', 'GeneID', 'TAIRLoc2',
                                     'Symbol1', 'Symbol2', 'NCBITaxo',
                                     'Source', 'Type', 'Annotation'))

Symbol1 <- athGORaw$Symbol1 %>%
  str_trim %>%
  str_split('\\|') %>%
  sapply(function(x) {
    x %<>%
      str_detect('AT\\dG') %>%
      x[.] %>%
      .[1]
    return(x)
  })

Symbol2 <- athGORaw$Symbol2 %>%
  if_else(str_detect(., 'AT\\dG\\d{5}') & !is.na(.), ., NA_character_)

Symbol <- sapply(seq_along(Symbol1), function(x) {

  eachSymbol <- c(Symbol1[x], Symbol2[x]) %>%
    {.[!is.na(.)][1]}

  return(eachSymbol)
})
################################################################

#########################GO from TAIR##############################
## ath GO is downloaded from TAIR at April 29, 2019.
library('readr')
library('dplyr')
library('stringr')

setwd('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/geneset')

gffAnno <- read_csv('../results/Ensembl_ath_Anno.csv',
                    col_types = cols(Chromosome = col_character()))

athGORaw <- read_delim('ATH_GO_GOSLIM.txt',
                       delim = '\t',
                       col_names = c('Symbol', 'Locus', 'TranscriptsID',
                                     'Anno1', 'Anno2', 'GO',
                                     'unknowNum', 'Aspect', 'Anno3',
                                     'Evidencewith', 'EvidenceDesc', 'Domains',
                                     'Ref', 'Source', 'Date'))

athGORaw %<>%
  select(Symbol, GO) %>%
  distinct %>%
  filter(Symbol %>% {str_detect(., 'AT(\\d|C|M)G\\d{5}') & !is.na(.)})

## upper, be careful at2g32273
athGO <- gffAnno %>%
  mutate(GeneID = ID %>%
           str_extract('(AT|at)(\\d|C|M)(G|g)\\d{5}') %>%
           str_to_upper) %>%
  select(ID, GeneID) %>%
  inner_join(athGORaw, by = c('GeneID' = 'Symbol')) %>%
  select(ID, GO) %>%
  {split(.$ID, .$GO)}

save(athGO, file = 'athGO.RData', compress = 'xz')
###################################################################


##########################KEGG pathway ################################
setwd('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/geneset/')

library('KEGGAPI')
library('readr')
library('tibble')
library('magrittr')
library('stringr')
library('dplyr')

gffAnno <- read_csv('../results/Ensembl_ath_Anno.csv',
                    col_types = cols(Chromosome = col_character()))

## get KEGG path
keggRaw <- getKEGGPathGenes('ath')
keggRaw <- sapply(keggRaw, function(x) {
  eachID <- sapply(strsplit(x, split = ':', fixed = TRUE), '[[', 2)
  return(eachID)
})
keggMat <- tibble(ID = unlist(keggRaw),
                  pathID = keggRaw %>% names %>% rep(sapply(keggRaw, length)))

## check later
## keggMat$ID %>%
##   str_detect('Arth') %>%
##   keggMat$ID[.] %>%
##   unique

athKEGG <- gffAnno %>%
  mutate(GeneID = ID %>%
           str_extract('(AT|at)(\\d|C|M)(G|g)\\d{5}') %>%
           str_to_upper) %>%
  select(ID, GeneID) %>%
  inner_join(keggMat, by = c('GeneID' = 'ID')) %>%
  select(ID, pathID) %>%
  {split(.$ID, .$pathID)}

save(athKEGG, file = 'athKEGG.RData', compress = 'xz')
#######################################################################


################################BioCyc genes#########################
setwd('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/geneset/')

library('BioCycAPI') ## version 0.2.1
library('ParaMisc')
library('doParallel')
library('foreach')
library('magrittr')

athcycFolder <- 'athcycgenes'

if (!dir.exists(athcycFolder)) {
  dir.create(athcycFolder)
} else {}

cycIDsRaw <- getCycGenes('ARA')

##~~~~~~~~~~~~parallel download~~~~~~~~~~~~~~~~~
registerDoParallel(cores = 12)

cutMat <- CutSeqEqu(length(cycIDsRaw), 12)

for (j in 1:ncol(cutMat)) {

  print(paste0('It is running ', j, ' in a total of ', ncol(cutMat), '.'))

  cycAnno <- foreach(i = cutMat[1, j] : cutMat[2, j]) %dopar% {
    eachCycIDs <- getCycGeneInfo(cycIDsRaw[i])
    return(eachCycIDs)
  }
  names(cycAnno) <- cycIDsRaw[cutMat[1, j] : cutMat[2, j]]

  paste0('ath', cutMat[1, j], '_', cutMat[2, j], '.RData') %>%
    file.path(athcycFolder, .) %>%
    save(cycAnno, file = ., compress = 'xz')
}

stopImplicitCluster()
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#######################################################################


#######################BioCyc pathways ath####################
setwd('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/geneset/')

library('BioCycAPI') ## version 0.2.1
library('ParaMisc')
library('doParallel')
library('foreach')
library('magrittr')

athcycFolder <- 'athcycpaths'

if (!dir.exists(athcycFolder)) {
  dir.create(athcycFolder)
} else {}

athPath <- getCycPathway('ARA')

##~~~~~~~~~~~~parallel download~~~~~~~~~~~~~~~~~
registerDoParallel(cores = 12)

cutMat <- CutSeqEqu(nrow(athPath), 12)

for (j in 1:ncol(cutMat)) {

  print(paste0('It is running ', j, ' in a total of ', ncol(cutMat), '.'))

  cycAnno <- foreach(i = cutMat[1, j] : cutMat[2, j]) %dopar% {
    eachCycIDs <- athPath[i, 1] %>%
      as.character %>%
      getCycGenesfPathway
    return(eachCycIDs)
  }
  names(cycAnno) <- athPath$pathID[cutMat[1, j] : cutMat[2, j]]

  paste0('ath', cutMat[1, j], '_', cutMat[2, j], '.RData') %>%
    file.path(athcycFolder, .) %>%
    save(cycAnno, file = ., compress = 'xz')
}

stopImplicitCluster()
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###########################################################

#########################merge ath biocyc path##################
setwd('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/geneset/')

library('KEGGAPI') ## version 0.1.7.4
library('BioCycAPI') ## version 0.2.1
library('readr')
library('magrittr')
library('foreach')
library('tibble')
library('dplyr')
library('stringr')

gffAnno <- read_csv('../results/Ensembl_ath_Anno.csv',
                    col_types = cols(Chromosome = col_character()))

## saocyc genes
athcycFolder <- 'athcycgenes'
athFiles <- dir(athcycFolder, full.names = TRUE)

athConv <- foreach(i = seq_along(athFiles), .combine = bind_rows) %do% {

  load(athFiles[i])

  eachConv <- tibble(
    cycID = names(cycAnno),
    ID1 = sapply(cycAnno, function(x) {
      l <- str_detect(x$name, '(AT|at)(\\d|C|M)(G|g)\\d{5}')
      eachname <- x$name[l][1] %>%
        str_extract('(AT|at)(\\d|C|M)(G|g)\\d{5}') %>%
        str_to_upper
      return(eachname)
    }))

  return(eachConv)
}

athConv %<>%
  mutate(ID2 = cycID %>% str_extract('(AT|at)(\\d|C|M)(G|g)\\d{5}'))

ID <- sapply(seq_along(athConv$ID1), function(x) {

  eachSymbol <- c(athConv$ID1[x], athConv$ID2[x]) %>%
    {.[!is.na(.)][1]}

  return(eachSymbol)
})

athConv %<>%
  mutate(ID = ID) %>%
  select(-ID1, -ID2) %>%
  filter(!is.na(ID))

## biocyc pathways
athcycFolder <- 'athcycpaths'
athFiles <- dir(athcycFolder, full.names = TRUE)

athPath <- foreach(i = seq_along(athFiles), .combine = append) %do% {
  load(athFiles[i])
  eachPath <- lapply(cycAnno, '[[', 1)
  return(eachPath)
}

pathMat <- tibble(ID = unlist(athPath),
                  pathID = athPath %>% names %>% rep(sapply(athPath, length))) %>%
  inner_join(athConv, by = c('ID' = 'cycID')) %>%
  select(-ID) %>%
  rename(ID = ID.y)

athBioCyc <- gffAnno %>%
  mutate(GeneID = ID %>%
           str_extract('(AT|at)(\\d|C|M)(G|g)\\d{5}') %>%
           str_to_upper) %>%
  select(ID, GeneID) %>%
  inner_join(pathMat, by = c('GeneID' = 'ID')) %>%
  select(ID, pathID) %>%
  {split(.$ID, .$pathID)}

save(athBioCyc, file = 'athBioCyc.RData', compress = 'xz')
################################################################
