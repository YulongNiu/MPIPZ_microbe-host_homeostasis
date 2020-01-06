####################merge fq files lines or runs##################
library('magrittr')
library('doParallel')
library('foreach')
library('readr')
library('dplyr')

rawfqPath <- '/biodata/dep_psl/grp_rgo/yniu/KaWai_soil_Syncom'
resFolder <- '/netscratch/dep_psl/grp_rgo/yniu/KaWaiFlg22/raw_data_soil'
catPath <- '/bin/cat'
mvPath <- '/bin/mv'
ncore <- 40

rawfq <- dir(rawfqPath,
             pattern = 'fastq.gz')

fqs <- rawfq %>%
  strsplit('_', fixed = TRUE) %>%
  lapply('[', c(1, 2, 7)) %>%
  sapply(paste, collapse = '_')

## group fq gz files
fqIdx <- split(seq_along(fqs), fqs)
fqPrefix <- names(fqIdx)

## check unique md5sum
registerDoParallel(cores = ncore)
fqmd5 <- foreach (i = seq_along(fqIdx)) %dopar% {
  eachmd5 <- fqIdx[[i]] %>%
    {file.path(rawfqPath, rawfq[.])} %>%
    paste('md5sum ', .) %>%
    system(intern = TRUE) %>%
    strsplit(split = ' ', fixed = TRUE) %>%
    unlist %>%
    .[c(1, 3)]

  return(eachmd5)
} %>% do.call(rbind, .)

(fqmd5[, 1] %>% unique %>% length) == nrow(fqmd5)
stopImplicitCluster()

## merge
registerDoParallel(cores = ncore)
foreach (i = seq_along(fqIdx), .combine = c) %dopar% {
  ## input files
  fqin <- fqIdx[[i]] %>%
    {file.path(rawfqPath, rawfq[.])} %>%
    paste(collapse = ' ')

  fqout <- fqPrefix[i] %>%
    paste0('.fq.gz') %>%
    file.path(resFolder, .)

  mergeC <- paste(catPath,
                  fqin,
                  '>',
                  fqout)

  print(mergeC)

  system(mergeC)

  return(NULL)
}
stopImplicitCluster()

##~~~~~~~~~~~~~~~~~~~~~raw~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## change to sample names
fqraws <- dir(resFolder)
fqnews <- c('Mock', 'Flg22', 'Flg22_SynCom33', 'Flg22_SynCom35') %>%
  rep(3) %>%
  paste0('_', rep(1:3, each = 4), '.fq.gz')

for (i in seq_along(fqraws)) {

  fqin <- fqraws[i] %>%
    file.path(resFolder, .)

  fqout <- fqnews[i] %>%
    file.path(resFolder, .)

  mvC <- paste(mvPath,
               fqin,
               fqout)

  print(mvC)

  system(mvC)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~1stadd~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
anno <- read_csv('/netscratch/dep_psl/grp_rgo/yniu/KaWaiFlg22/results/Ka-Wai_sample_4313.csv')
fqraws <- dir(resFolder)

for (i in seq_len(nrow(anno))) {

  ## merge different batch
  fqin <- anno[i, c(1, 6)] %>%
    as.character %>%
    .[!is.na(.)] %>%
    file.path(resFolder, .) %>%
    paste0('_R1.fq.gz') %>%
    paste(collapse = ' ')

  fqout <- anno[i, 5] %>%
    file.path(resFolder, .) %>%
    paste0('_R1.fq.gz')

  mergeC <- paste(catPath,
                  fqin,
                  '>',
                  fqout)
  print(mergeC)

  system(mergeC)

  fqin <- anno[i, c(1, 6)] %>%
    as.character %>%
    .[!is.na(.)] %>%
    file.path(resFolder, .) %>%
    paste0('_R2.fq.gz') %>%
    paste(collapse = ' ')

  fqout <- anno[i, 5] %>%
    file.path(resFolder, .) %>%
    paste0('_R2.fq.gz')

  mergeC <- paste(catPath,
                  fqin,
                  '>',
                  fqout)
  print(mergeC)

  system(mergeC)
}

system('rm `ls | grep "4313\\|4219"`')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~soil~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
anno <- read_csv('/netscratch/dep_psl/grp_rgo/yniu/KaWaiFlg22/results/Ka-Wai_soil.csv')
fqraws <- dir(resFolder)

for (i in seq_len(nrow(anno))) {

  ## merge different batch
  fqin <- anno[i, c(1, 6)] %>%
    as.character %>%
    .[!is.na(.)] %>%
    file.path(resFolder, .) %>%
    paste0('_R1.fq.gz') %>%
    paste(collapse = ' ')

  fqout <- anno[i, 5] %>%
    file.path(resFolder, .) %>%
    paste0('_R1.fq.gz')

  mergeC <- paste(catPath,
                  fqin,
                  '>',
                  fqout)
  print(mergeC)

  system(mergeC)

  fqin <- anno[i, c(1, 6)] %>%
    as.character %>%
    .[!is.na(.)] %>%
    file.path(resFolder, .) %>%
    paste0('_R2.fq.gz') %>%
    paste(collapse = ' ')

  fqout <- anno[i, 5] %>%
    file.path(resFolder, .) %>%
    paste0('_R2.fq.gz')

  mergeC <- paste(catPath,
                  fqin,
                  '>',
                  fqout)
  print(mergeC)

  system(mergeC)
}

system('rm `ls | grep "4296\\|4397"`')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################################################################
