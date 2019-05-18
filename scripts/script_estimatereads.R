#####################estimate from transcripts######################
library('stringr')
library('utils')
library('readr')
library('dplyr')
library('magrittr')
library('tibble')

gffPath <- '/extDisk1/Biotools/RefData/ath/Arabidopsis_thaliana.TAIR10.43.gff3.gz'
gffAnno <- read_tsv(gffPath,
                    col_names = c('chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'),
                    col_types = cols(chromosome = col_character()),
                    comment = '#')

##~~~~~~~~~~exon~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
exonAnno <- gffAnno %>%
  filter(feature %in% c('exon')) %>%
  select(chromosome, source, feature, start, end) %>%
  distinct %>% ## remove repeat exons
  mutate(length = abs(end - start + 1))

readlen <- 150

## depth
depth <- c(0.1, 0.5, 1, 2, 3, 5, 10, 15, 20, 30, 40, 50)

tibble(Coverage = paste0(depth, 'X'),
       `Reads (million)` = sum(exonAnno$length) / 1e6 / readlen * depth)

tibble(Coverage = paste0(depth, 'X'),
       `Reads (million)` = 27655 * 6.7 * 335.5 / 1e6 / readlen * depth)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~transcripts~~~~~~~~~~~~~~~~~~~~~~~
geneAnno <- gffAnno %>%
  filter(feature %in% c('lnc_RNA', 'miRNA', 'mRNA', 'ncRNA',
                        'rRNA', 'snoRNA', 'snRNA', 'tRNA'))

noteAnno <- geneAnno %>%
  select(attribute) %>%
  unlist %>%
  unname %>%
  str_trim

geneAnno %<>%
  mutate(BioType = str_extract(noteAnno, '(?<=biotype=).*?(?=;)')) %>%
  filter(BioType %in% c('protein_coding', 'lncRNA')) %>%
  select(chromosome, source, feature, start, end) %>%
  distinct %>% ## remove repeat transcripts
  mutate(length = abs(end - start + 1))

readlen <- 150

## depth
depth <- c(0.1, 0.5, 1, 2, 3, 5, 10, 15, 20, 30, 40, 50)

tibble(Coverage = paste0(depth, 'X'),
       `Reads (million)` = sum(geneAnno$length) / 1e6 / readlen * depth)



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################################################################
