###########################Raw reads######################
rawpath <- '/netscratch/dep_psl/grp_rgo/yniu/KaWaiFlg22/raw_data'
setwd(rawpath)

library('magrittr')
library('doParallel')
library('foreach')
library('tibble')
library('readr')

ncore <- 12

fqs <- dir(rawpath,
           pattern = 'fq.gz',
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
  sapply('[[', 1)

tibble(sample = snames,
       rawfq = rn) %>%
  write_csv('raw_seqnumber.csv')
##########################################################

#################extract Kallisto and HISAT2 output########
setwd('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/')

library('dplyr')
library('readr')

KHoutput <- function(op, type = 'PE', org = 'hsa') {

  ## INPUT: 'op' is a character vector. 'type' is 'PE' (pair-end) or 'SE' (single-end). 'org' is the organism name.
  ## OUTPUT: A tibble, 1st column is input reads number, 2nd column is Kallisto aligned read number, and 3rd column is the HISAT2 aligned read number.
  ## USAGE: Extract the number of aligned reads.

  require('stringr')
  require('magrittr')
  require('tibble')

  ## input reads number
  fqnum <- op %>%
    str_detect('\\d+ reads; of these:') %>%
    op[.] %>%
    strsplit(split = ' ', fixed = TRUE) %>%
    sapply('[[', 1) %>%
    as.numeric

  ## HISAT2 aligned
  hmapnum <- op %>%
    str_detect('.* aligned 0 times$') %>%
    op[.] %>%
    strsplit(split = ' ', fixed = TRUE) %>%
    {  if (type == 'PE') {
         sapply(., '[[', 9)
       } else {
         sapply(., '[[', 5)
       }} %>%
    as.numeric %>%
    {  if (type == 'PE') {
         ./2
       } else .} %>%
    {fqnum - .}

  ## Kallisto aligned
  kmapnum <- op %>%
    str_detect('.* reads pseudoaligned') %>%
    op[.] %>%
    strsplit(split = ' ', fixed = TRUE) %>%
    sapply('[[', 5) %>%
    str_replace_all(',', '') %>%
    as.numeric

  ## sample names
  snames <- op %>%
    str_detect('HISAT2 using') %>%
    op[.] %>%
    strsplit(split = ' ', fixed = TRUE) %>%
    sapply('[[', 6) %>%
    {substr(., 1, nchar(.) - 1)}

  res <- tibble(sample = snames,
                trimfq = fqnum,
                hmap = hmapnum,
                kmap = kmapnum,
                org = org)

  return(res)
}

##~~~~~~~~~~~~~~~~~~~~~~~test contamination~~~~~~~~~~~~~~~~~~~~~~~~~~~
athout <- 'align_nohup.out' %>%
  readLines %>%
  KHoutput(type = 'SE', org = 'ath') %>%
  mutate(H_ath = round(hmap/trimfq, 3), K_ath = round(kmap/trimfq, 3)) %>%
  select(c(-hmap, -kmap, -org))

## raw reads
rawrd <- read_csv('raw_seqnumber.csv')

contam <- rawrd %>%
  inner_join(athout) %>%
  slice(c(10:12, 1:9))

write_csv(contam, 'ath_alignment.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
######################################################
