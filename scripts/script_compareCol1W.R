##########################Compare Col1 and Wmutant#####################
library('magrittr')
library('readr')
library('dplyr')
library('VennDiagram')

setwd('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results')

col1Raw <- read_csv('cluster10_g4_1stadd.csv',
                    col_types = cols(Chromosome = col_character()))
wmutantRaw <- read_csv('kmeans_10.csv',
                       col_types = cols(Chromosome = col_character()))

cgenes <- c('AT1G14550.1', 'AT2G30750.1', 'AT2G19190.1', 'AT5G24110.1')

col1Raw %>%
  filter(ID %in% cgenes) %>%
  .$cl

wmutantRaw %>%
  filter(ID %in% cgenes) %>%
  .$cl

## wmutant cluster 9
## col1 cluster 10

col110ID <- col1Raw %>%
  filter(cl == 10) %>%
  .$ID

wmutant9ID <- wmutantRaw %>%
  filter(cl == 9) %>%
  .$ID

cairo_pdf('veen_plot.pdf')
grid.newpage()
draw.pairwise.venn(setdiff(col110ID, wmutant9ID) %>% length,
                   setdiff(wmutant9ID, col110ID) %>% length,
                   intersect(col110ID, wmutant9ID) %>% length,
                   category = c('Col0', 'Mutant'),
                   lty = rep('blank', 2),
                   fill = c('light blue', 'pink'),
                   alpha = rep(0.5, 2),
                   cat.pos = c(0, 0),
                   cat.dist = rep(0.025, 2))
dev.off()

intersect(col110ID, wmutant9ID) %>%
  {filter(wmutantRaw, ID %in% .)} %>%
  mutate_at(c('Gene', 'Description'), .funs = list(~if_else(is.na(.), '', .))) %>%
  write_csv('tmp1.csv')

########################################################################
