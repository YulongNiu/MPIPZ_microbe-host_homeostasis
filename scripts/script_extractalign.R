###########################Raw reads######################
library('magrittr')
library('doParallel')
library('foreach')
library('tibble')
library('readr')

ncore <- 25

rawpath <- '/netscratch/dep_psl/grp_rgo/yniu/AmeliaMaize/raw_data'

fqs <- dir(rawpath,
           pattern = 'R1.fq.gz',
           full.names = TRUE)

registerDoParallel(cores = ncore)

rn <- foreach(i = seq_along(fqs), .combine = c) %dopar% {

  eachrn <- paste('zcat', fqs[i], '| awk "END{print NR/4}"') %>%
    system(inter = TRUE) %>%
    as.numeric

  return(eachrn)
}

stopImplicitCluster()

snames <- fqs %>%
  strsplit(split = '/', fixed = TRUE) %>%
  sapply('[[', 8) %>%
  strsplit(split = '.', fixed = TRUE) %>%
  sapply('[[', 1) %>%
  {substr(., 1, nchar(.) - 3)}

tibble(sample = snames,
       rawfq = rn) %>%
  write_csv('raw_seqnumber.csv')
##########################################################

###########################################################
