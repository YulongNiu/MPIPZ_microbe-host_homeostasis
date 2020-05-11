
## originally by Yulong Niu
## yulong.niu@hotmail.com

##########################Compare Col0 and Wmutant#####################
library('magrittr')
library('readr')
library('dplyr')
library('VennDiagram')

setwd('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results')

col0Raw <- read_csv('cluster10_g4_1stadd.csv',
                    col_types = cols(Chromosome = col_character()))
wmutantRaw <- read_csv('kmeans_10.csv',
                       col_types = cols(Chromosome = col_character()))

cgenes <- c('AT1G14550.1', 'AT2G30750.1', 'AT2G19190.1', 'AT5G24110.1')

col0Raw %>%
  filter(ID %in% cgenes) %>%
  .$cl

wmutantRaw %>%
  filter(ID %in% cgenes) %>%
  .$cl

##    -
## _ /  \_
## wmutant cluster 9
## col0 cluster 10

## _   /
##  \ /
##   -
## wmutant cluster 4
## col0 cluster 1

col010ID <- col0Raw %>%
  filter(cl == 10) %>%
  .$ID

wmutant9ID <- wmutantRaw %>%
  filter(cl == 9) %>%
  .$ID

cairo_pdf('veen_plot_35up.pdf')
grid.newpage()
draw.pairwise.venn(setdiff(col010ID, wmutant9ID) %>% length,
                   setdiff(wmutant9ID, col010ID) %>% length,
                   intersect(col010ID, wmutant9ID) %>% length,
                   category = c('Col0', 'Mutant'),
                   lty = rep('blank', 2),
                   fill = c('light blue', 'pink'),
                   alpha = rep(0.5, 2),
                   cat.pos = c(0, 0),
                   cat.dist = rep(0.025, 2))
dev.off()

intersect(col010ID, wmutant9ID) %>%
  {filter(wmutantRaw, ID %in% .)} %>%
  mutate_at(c('Gene', 'Description'), .funs = list(~if_else(is.na(.), '', .))) %>%
  write_csv('inter_col0mu_35down.csv')

intersect(col010ID, wmutant9ID) %>%
  {filter(col0Raw, ID %in% .)} %>%
  mutate_at(c('Gene', 'Description'), .funs = list(~if_else(is.na(.), '', .))) %>%
  write_csv('inter_col0mu_35down_1stadd.csv')
########################################################################
