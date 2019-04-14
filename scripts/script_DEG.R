###########################DEGs##################################
##~~~~~~~~~~~~~~~~~~~~load k alignments~~~~~~~~~~~~~~~~~~~~~~~~~~
library('tximport')
library('rhdf5')
library('magrittr')

wd <- '/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/bamfiles'
setwd(wd)


athlabel <- rep(c('Mock', 'Flg22', 'Flg22_SynCom33', 'Flg22_SynCom35'), each = 3) %>%
  paste(rep(1:3, 4), sep = '_') %>%
  paste0('_ath_kallisto')

files <- file.path(wd, athlabel, 'abundance.h5')
names(files) <- rep(c('Mock', 'Flg22', 'Flg22_SynCom33', 'Flg22_SynCom35'), each = 3) %>%
  paste(rep(1:3, 4), sep = '_')
athk <- tximport(files, type = 'kallisto', txOut = TRUE)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~DEG anlaysis~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################
