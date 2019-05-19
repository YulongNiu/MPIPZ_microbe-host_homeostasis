#####################estimate from transcripts######################
library('stringr')
library('utils')
library('readr')
library('dplyr')
library('magrittr')
library('tibble')
library('tximport')
library('rhdf5')
library('ggplot2')


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

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~real data~~~~~~~~~~~~~~~~~~~~~~~~~
wd <- '/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/bamfiles'
setwd(wd)

slabel <- rep(c('Mock', 'Flg22', 'Flg22_SynCom33', 'Flg22_SynCom35'), each = 3) %>%
  paste(rep(1:3, 4), sep = '_') %>%
  paste0('_ath_kallisto')

files <- file.path(wd, slabel, 'abundance.h5')
names(files) <- rep(c('Mock', 'Flg22', 'Flg22_SynCom33', 'Flg22_SynCom35'), each = 3) %>%
  paste(rep(1:3, 4), sep = '_')
kres <- tximport(files, type = 'kallisto', txOut = TRUE)

kres$counts %>%
  as.data.frame %>%
  as_tibble %>%
  filter(Mock_1 > 0) %>%
  ggplot(aes(.$Mock_1)) +
  geom_histogram(aes(y =..density..))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################################################################

#######################subSeq#######################################
library('subSeq')
library('magrittr')
library('DESeq2')

setwd('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/')

load('degres_condi_Mock.RData')

rawCount <- counts(degres)

tvscIdx <- c(1:6, c(1:3, 7:9), c(1:3, 10:12))
tvscName <- c('Flg22_vs_Mock_subSeq', 'Flg22_SynCom33_vs_Mock_subSeq', 'Flg22_SynCom35_vs_Mock_subSeq')

for(i in seq_along(tvscIdx)) {
  ## two conditions
  testCounts <- rawCount[, 1:6] %>%
    .[rowSums(.) >= 5, ]
  testCondi <- degres$condition %>%
    .[1:6] %>%
    droplevels

  proportions <- 10^seq(-2, 0, .1)
  ss <- subsample(testCounts, proportions, method = c('edgeR', 'voomLimma', 'DESeq2'), treatment = testCondi)

  seed <- getSeed(ss)
  ssglm <- subsample(testCounts, proportions, method = c('edgeR.glm'), mm= model.matrix(~ testCondi), seed = seed)
  ssComb <- combineSubsamples(ss, ssglm)

  plot(ssComb)
  ggsave(paste0(tvscName[i], '.pdf'))
  ggsave(paste0(tvscName[i], '.jpg'))
}


###################################################################
