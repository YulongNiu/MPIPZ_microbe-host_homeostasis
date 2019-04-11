####################merge fq files lines or runs##################
require('magrittr')
require('doParallel')
require('foreach')

rawfqPath <- '/biodata/dep_psl/grp_rgo/metatranscriptomics/data/flg22'
resFolder <- '/netscratch/dep_psl/grp_rgo/yniu/KaWaiFlg22/raw_data'
zcatPath <- '/bin/cat'
ncore <- 12

rawfq <- dir(rawfqPath,
             pattern = 'fastq.gz')

fqs <- rawfq %>%
  strsplit('_', fixed = TRUE) %>%
  lapply('[', c(1, 2, 4, 7)) %>%
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

  mergeC <- paste(zcatPath,
                  fqin,
                  '>',
                  fqout)

  print(mergeC)

  system(mergeC)

  return(NULL)
}

stopImplicitCluster()


####################################################################
