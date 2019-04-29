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

setwd('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/')

gffAnno <- read_csv('results/Ensembl_ath_Anno.csv',
                    col_types = cols(Chromosome = col_character()))

athGORaw <- read_delim('ATH_GO_GOSLIM.txt',
                       delim = '\t',
                       col_names = c('Symbol', 'Locus', 'TranscriptsID',
                                     'Anno1', 'Anno2', 'GO',
                                     'unknowNum', 'Aspect', 'Anno3',
                                     'Evidencewith', 'EvidenceDesc', 'Domains',
                                     'Ref', 'Source', 'Date'))


athGO <- athGORaw %>%
  select(Symbol, GO) %>%
  distinct %>%
  filter(Symbol %>% {str_detect(., 'AT\\dG\\d{5}') & !is.na(.)})

athGO <- gffAnno %>%
  mutate(GeneID = ID %>% str_extract('AT\\dG\\d{5}')) %>%
  select(ID, GeneID) %>%
  inner_join(athGO, by = c('GeneID' = 'Symbol')) %>%
  select(ID, GO) %>%
  {split(.$ID, .$GO)}

save(athGO, file = 'geneset/athGO.RData', compress = 'xz')
###################################################################
