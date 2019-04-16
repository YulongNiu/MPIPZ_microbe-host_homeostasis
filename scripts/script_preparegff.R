########################deal with gff3########################
library('stringr')
library('utils')
library('readr')
library('dplyr')
library('magrittr')

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
                    Gene = str_extract(noteAnno, '(?<=Name=).*?(?=;)') %>% {if_else(is.na(.), '', .)} %>% sapply(URLdecode),
                    BioType = str_extract(noteAnno, '(?<=biotype=).*?(?=;)'),
                    Description = str_extract(noteAnno, '(?<=description=).*?(?=;)') %>% {if_else(is.na(.), '', .)} %>% sapply(URLdecode))
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~cDNA table~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cDNAanno <- gffAnno %>%
  filter(feature %in% c('lnc_RNA', 'miRNA', 'mRNA', 'ncRNA',
                        'rRNA', 'snoRNA', 'snRNA', 'tRNA'))
## "ID=transcript:ATCG00310.1;Parent=gene:ATCG00310;biotype=tRNA;transcript_id=ATCG00310.1"
noteAnno <- cDNAanno %>%
  select(attribute) %>%
  unlist %>%
  unname %>%
  str_trim

cDNATable <- tibble(cDNA = str_extract(noteAnno, '(?<=ID=transcript:).*?(?=;)'),
                    ID = str_extract(noteAnno, '(?<=Parent=gene:).*?(?=;)'),
                    Chromosome = cDNAanno$chromosome,
                    Start = cDNAanno$start,
                    End = cDNAanno$end,
                    Strand = cDNAanno$strand)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## merge cDNA table and gene table
athAnno <- inner_join(geneTable, cDNATable, by = 'ID') %>%
  select(-ID) %>%
  select(cDNA, everything()) %>%
  rename(ID = cDNA) %>%
  mutate(Length = abs(End - Start)) %>%
  select(ID, Gene, Chromosome:Length, BioType, Description)

write_csv(athAnno, '/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/Ensembl_ath_Anno.csv')

##~~~~~~~~~~~~~~~~~~check cDNA in k~~~~~~~~~~~~~~~~~~~~~~~~~~~
kcDNA <- read_tsv('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/bamfiles/Flg22_1_ath_kallisto/abundance.tsv')

## kcDNA 48359
## athAnno 54013

## all kcDNA in athAnno
sum(kcDNA$target_id %in% athAnno$cDNA)

## lncRNA  miRNA  ncRNA   rRNA snoRNA  snRNA   tRNA
##   3879    325    377     15    287     82    689
anti_join(athAnno, kcDNA, by = c('cDNA' = 'target_id')) %>%
  .$BioType %>%
  table
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##############################################################
