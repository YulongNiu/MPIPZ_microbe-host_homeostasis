####################merge fq files lines or runs##################
library('magrittr')
library('doParallel')
library('foreach')
library('readr')
library('dplyr')

rawfqPath <- '/biodata/dep_psl/grp_rgo/yniu/KaWai_raw_data_1stadd'
resFolder <- '/netscratch/dep_psl/grp_rgo/yniu/KaWaiFlg22/raw_data_1stadd_test'
catPath <- '/bin/cat'
mvPath <- '/bin/mv'
ncore <- 12

rawfq <- dir(rawfqPath,
             pattern = 'fastq.gz')

fqs <- rawfq %>%
  strsplit('_', fixed = TRUE) %>%
  lapply('[', c(1, 2, 7)) %>%
  sapply(paste, collapse = '_')

## group fq gz files
fqIdx <- split(seq_along(fqs), fqs)
fqPrefix <- names(fqIdx)

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
anno <- read_delim('/netscratch/dep_psl/grp_rgo/yniu/KaWaiFlg22/results/list_samples_1stadd.txt', delim = '\t')
fqraws <- dir(resFolder)

for (i in seq_len(nrow(anno))) {

  fqin <- anno[i, 1] %>%
    file.path(resFolder, .) %>%
    paste0('_R1.fq.gz')

  fqout <- anno[i, 5] %>%
    file.path(resFolder, .) %>%
    paste0('_R1.fq.gz')

  mvC <- paste(mvPath,
               fqin,
               fqout)

  print(mvC)

  system(mvC)

  fqin <- anno[i, 1] %>%
    file.path(resFolder, .) %>%
    paste0('_R2.fq.gz')

  fqout <- anno[i, 5] %>%
    file.path(resFolder, .) %>%
    paste0('_R2.fq.gz')

  mvC <- paste(mvPath,
               fqin,
               fqout)

  print(mvC)

  system(mvC)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################################################################
