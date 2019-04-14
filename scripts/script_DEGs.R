########################deal with gff3########################
library('stringr')
library('utils')
library('foreach')
library('doParallel')
library('readr')
library('dplyr')
library('magrittr')

ncore <- 12

registerDoParallel(cores = ncore)

gffPath <- '/extDisk1/Biorefs/Arabidopsis_thaliana.TAIR10.42.gff3.gz'
gffAnno <- read_tsv(gffPath,
                    col_names = c('chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'),
                    col_types = cols(chromosome = col_character()),
                    comment = '#')

##~~~~~~~~~~~~~~~~~~~~~~~~~~gene and lncgene table~~~~~~~~~~~~~~
geneAnno <- gffAnno %>%
  filter(feature %in% c('gene', 'ncRNA_gene'))

## "ID=gene:AT1G01760;biotype=protein_coding;description=TAD1 [Source:UniProtKB/TrEMBL%3BAcc:A0A178W782];gene_id=AT1G01760;logic_name=araport11"
noteAnno <- geneAnno %>%
  select(attribute) %>%
  unlist %>%
  unname %>%
  str_trim

geneTable <- tibble(ID = str_extract(noteAnno, '(?<=ID=gene:).*?(?=;)'),
                    Gene = str_extract(noteAnno, '(?<=Name=).*?(?=;)') %>% {if_else(is.na(.), '', .)},
                    Chromosome = geneAnno$chromosome,
                    Start = geneAnno$start,
                    End = geneAnno$end,
                    Strand = geneAnno$strand,
                    BioType = str_extract(noteAnno, '(?<=biotype=).*?(?=;)'),
                    Description = str_extract(noteAnno, '(?<=description=).*?(?=;)') %>% {if_else(is.na(.), '', .)} %>% sapply(URLdecode))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




%>%
  filter(feature %in% c('lnc_RNA', 'miRNA', 'mRNA', 'rRNA', 'snoRNA', 'snRNA'))

noteAnno <- gffAnno %>%
  select(attribute) %>%
  unlist %>%
  unname %>%
  str_trim

ids <- noteAnno %>%
  strsplit(., split = ';', fixed = TRUE) %>%
  sapply(., '[[', 1) %>%
  strsplit(., split = ':', fixed = TRUE) %>%
  sapply(., '[[', 2)

##"ID=gene:AT1G01760;biotype=protein_coding;description=TAD1 [Source:UniProtKB/TrEMBL%3BAcc:A0A178W782];gene_id=AT1G01760;logic_name=araport11"
ExtractNote <- function(x) {
  x <- str_extract(x, 'Note=.*')
  ## check if no "Note=.*"
  if (is.na(x)) {
    return(x)
  } else {}

  x <- substring(x, 6)

  ## split with ';'
  xList <- unlist(str_split(x, ';'))
  eqIdx <- str_detect(xList, '=')
  if (sum(eqIdx) >= 1) {
    eq1stNum <- which(eqIdx)[1]
    x <- paste(xList[1:(eq1stNum-1)], collapse = ';')
  } else {}

  return(x)
}

ExtractGeneName <- function(x) {
  x <- str_extract(x, 'Gene=.*?;');

  if (!is.na(x)) {
    xLen <- nchar(x);
    x <- substr(x, 6, xLen - 1);
  } else {}

  return(x)
}

geneNames <- foreach(i = 1:length(noteAnno), .combine = c) %dopar% {
  x <- URLdecode(noteAnno[i])
  x <- ExtractGeneName(x)

  if (is.na(x)) {
    x <- ids[i]
  } else {}

  return(x)
}

noteAnno <- foreach(i = 1:length(noteAnno), .combine = c) %dopar% {
  x <- URLdecode(noteAnno[i])
  x <- ExtractNote(x)
  return(x)
}

gffAnno %<>%
  select(-attribute) %>%
  mutate(ids = ids, geneNames = geneNames, noteAnno = noteAnno)


kres <- read_tsv('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/bamfiles/Flg22_1_ath_kallisto/abundance.tsv')

##############################################################
