
## originally by Yulong Niu
## yulong.niu@hotmail.com

###########################Raw reads######################
rawpath <- '/netscratch/dep_psl/grp_rgo/yniu/KaWaiFlg22/raw_data_soil'
setwd(rawpath)

library('magrittr')
library('doParallel')
library('foreach')
library('tibble')
library('readr')
library('dplyr')

ncore <- 40

fqs <- dir(rawpath,
           pattern = 'fq.gz',
           full.names = TRUE)

registerDoParallel(cores = ncore)

rn <- foreach(i = seq(1, length(fqs), 2), .combine = bind_rows) %dopar% {

  eachr1 <- paste('zcat', fqs[i], '| awk "END{print NR/4}"') %>%
    system(inter = TRUE) %>%
    as.numeric

  eachr2 <- paste('zcat', fqs[i + 1], '| awk "END{print NR/4}"') %>%
    system(inter = TRUE) %>%
    as.numeric

  eachrn <- tibble(rawfq_R1 = eachr1, rawfq_R2 = eachr2)

  return(eachrn)
}

stopImplicitCluster()

snames <- fqs %>%
  strsplit(split = '/', fixed = TRUE) %>%
  sapply('[[', 8) %>%
  strsplit(split = '.', fixed = TRUE) %>%
  sapply('[[', 1) %>%
  substr(start = 1, stop = nchar(.) - 3) %>%
  unique

rn %>%
  mutate(sample = snames) %>%
  write_csv('../results/raw_seqnumber_soil.csv')
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
athout <- 'align_nohup_soil.out' %>%
  readLines %>%
  KHoutput(type = 'PE', org = 'ath') %>%
  mutate(H_ath = round(hmap/trimfq, 3), K_ath = round(kmap/trimfq, 3)) %>%
  select(c(-hmap, -kmap, -org))

athvirusout <- 'align_nohup_soil_latenvirus.out' %>%
  readLines %>%
  KHoutput(type = 'PE', org = 'ath') %>%
  mutate(H_ath = round(hmap/trimfq, 3), K_ath = round(kmap/trimfq, 3)) %>%
  select(c(-hmap, -kmap, -org)) %>%
  rename(Hvirus_ath = H_ath, Kvirus_ath = K_ath)

## raw reads
rawrd <- read_csv('raw_seqnumber_soil.csv')

contam <- rawrd %>%
  inner_join(athout) %>%
  inner_join(athvirusout) %>%
  select(sample, everything())

write_csv(contam, 'ath_alignment_soil.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
######################################################
