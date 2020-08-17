## originally by Yulong Niu
## yulong.niu@hotmail.com

###########################DEGs##################################
##~~~~~~~~~~~~~~~~~~~~~~~~~~~useful funcs~~~~~~~~~~~~~~~~~~~~~~~~
checkZeros <- function(v, threshold) {
  res <- ifelse(sum(v == 0) > threshold, FALSE, TRUE)
  return(res)
}

checkFlg22 <- function(v, threshold) {

  require('magrittr')

  res <- v %>%
    split(rep(1 : 10, each = 4)) %>%
    sapply(checkZeros, threshold) %>%
    all

  return(res)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


library('tximport')
library('rhdf5')
library('magrittr')
library('DESeq2')
library('tidyverse')
library('apeglm')

anno <- read_csv('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/Ensembl_ath_Anno.csv',
                 col_types = cols(Chromosome = col_character())) %>%
  mutate(Gene = Gene %>% {if_else(is.na(.), '', .)}) %>%
  mutate(Description = Description %>% {if_else(is.na(.), '', .)})


##~~~~~~~~~~~~~~~~~~~~load k alignments~~~~~~~~~~~~~~~~~~~~~~~~~~
wd <- '/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/flg22_crossref/Stringlis_2018/align_data'
setwd(wd)

sampleAnno <- c('Mock', 'flg22Pa', 'flg22417', 'WCS417', 'chitin') %>%
  rep(each = 12) %>%
  paste(rep(c('0p5h', '1h', '3h', '6h'), each = 3), sep = '_') %>%
  paste(rep(1:3, 16), sep = '_')
slabel <- sampleAnno %>%
  paste0('_kallisto')

files <- file.path(wd, slabel, 'abundance.h5')
names(files) <- sampleAnno
kres <- tximport(files, type = 'kallisto', txOut = TRUE)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~normalization~~~~~~~~~~~~~~~~~~~~~~~
setwd('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/flg22_crossref/Stringlis_2018/results/')

## sampleTable
condi <-  c('Mock', 'flg22Pa', 'flg22417', 'WCS417', 'chitin') %>%
  rep(each = 12) %>%
  paste(rep(c('0p5h', '1h', '3h', '6h'), each = 3), sep = '_')
sampleTable <- data.frame(condition = factor(condi, levels = condi %>% unique))
rownames(sampleTable) <- colnames(kres$counts)

degres <- DESeqDataSetFromTximport(kres, sampleTable, ~condition)

## remove 0|0|x, 0|0|0
degres %<>%
  estimateSizeFactors %>%
  counts(normalized = TRUE) %>%
  apply(1, checkFlg22, 1) %>%
  degres[., ]

degres <- DESeq(degres)

## count transformation
rld <- rlog(degres)
ntd <- normTransform(degres)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ##~~~~~~~~~~~~~~~~~~~~~~~~~~hidden batch effect~~~~~~~~~~~~~~~~~~~~~
## library('sva')
## library('ggplot2')

## dat <- rld %>%
##   assay %>%
##   {.[rowMeans(.) > 1, ]}
## mod <- model.matrix(~ condition, colData(degres))
## mod0 <- model.matrix(~ 1, colData(degres))

## ## manual detect surrogate variance
## svnum <- 4
## svseq <- svaseq(dat, mod, mod0, n.sv = svnum)

## ## auto detect sv
## svobj <- sva(dat, mod, mod0)
## svnum <- svobj$sv %>% ncol

## svobj$sv %>%
##   set_colnames(paste0('sv', seq_len(svnum))) %>%
##   as_tibble %>%
##   gather(key = 'sv', value = 'value') %>%
##   mutate(condition = colData(degres) %>%
##            .$condition %>%
##            rep(svnum) %>%
##            as.character,
##          sample = rep(colnames(degres), svnum)) %>%
##   mutate(group = paste(sv, condition, sep = '_')) %>%
##   ggplot(aes(sample, value, colour = sv, group = group)) +
##   geom_point() +
##   geom_line() +
##   theme(axis.text.x = element_text(angle = 90))
## ggsave('auto_sv_1stadd.jpg')
## ggsave('auto_sv_1stadd.pdf')
## ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~DEG Mock vs. 3 conditions~~~~~~~~~~~~~~~~~
## degres$sv1 <- svobj$sv[, 1]
## degres$sv2 <- svobj$sv[, 2]
## degres$sv3 <- svobj$sv[, 3]
## degres$sv4 <- svobj$sv[, 4]
## design(degres) <- ~sv1 + sv2 + sv3 + sv4 + condition

## degres <- DESeq(degres)

## each time point comparison
cond1 <- list(c('flg22Pa_0p5h', 'Mock_0p5h'),
              c('flg22Pa_1h', 'Mock_1h'),
              c('flg22Pa_3h', 'Mock_3h'),
              c('flg22Pa_6h', 'Mock_6h'),
              c('flg22417_0p5h', 'Mock_0p5h'),
              c('flg22417_1h', 'Mock_1h'),
              c('flg22417_3h', 'Mock_3h'),
              c('flg22417_6h', 'Mock_6h'),
              c('WCS417_0p5h', 'Mock_0p5h'),
              c('WCS417_1h', 'Mock_1h'),
              c('WCS417_3h', 'Mock_3h'),
              c('WCS417_6h', 'Mock_6h'))

cond <- c(cond1) %>%
  unique

resRaw <- lapply(cond,
                 function(x) {
                   degres %>%
                     results(contrast = c('condition', x)) %T>%
                     ## lfcShrink(coef = paste0('condition_', x), type = 'apeglm') %T>%
                     summary %>%
                     as_tibble %>%
                     select(pvalue, padj, log2FoldChange) %>%
                     rename_all(.funs = list(~paste0(paste(x, collapse = '_vs_'), '_', .)))
                 }) %>%
  bind_cols

res <- cbind.data.frame(as.matrix(mcols(degres)[, 1:10]), counts(degres, normalize = TRUE), stringsAsFactors = FALSE) %>%
  rownames_to_column(., var = 'ID') %>%
  as_tibble %>%
  bind_cols(resRaw) %>%
  inner_join(anno, by = 'ID') %>%
  select(ID, Gene : Description,  Mock_0p5h_1 : WCS417_6h_vs_Mock_6h_log2FoldChange) %>%
  arrange(flg22Pa_0p5h_vs_Mock_0p5h_padj)

write_csv(res, 'eachGroup_DEGs_Stringlis.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PCA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library('ggrepel')
library('ggplot2')
library('RColorBrewer')
library('limma')
library('sva')

## remove low count
dat <- rld %>%
  assay %>%
  {.[rowMeans(.) > 1, ]}

## group <- sampleTable$condition
## design <- model.matrix(~ group)
## rldData <- dat %>%
##   removeBatchEffect(covariates = svobj$sv,
##                     design = design)
rldData <- dat

## 1 - 2 C
pca <- prcomp(t(rldData))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
percentVar <- round(100 * percentVar)
pca1 <- pca$x[,1]
pca2 <- pca$x[,2]
pcaData <- tibble(PC1 = pca1, PC2 = pca2, Group = colData(rld)[, 1], ID = rownames(colData(rld))) %>%
  mutate(Treatment = rep(c('Mock', 'Mock+flg22Pa', 'Mock+flg22417', 'Mock+WCS417', 'Mock+chitin'), each = 12) %>% factor) %>%
  mutate(Time = rep(c('0.5h', '1h', '3h', '6h'), each = 3) %>% rep(5) %>% factor)

ggplot(pcaData, aes(x = PC1, y = PC2, colour = Treatment)) +
  geom_point(aes(shape = Time), size = 4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_colour_manual(values = brewer.pal(5, name = 'Set1')) +
  scale_shape_manual(values = c(17, 18, 15, 16),
                     name = 'Time point') +
  stat_ellipse(aes(x = PC1, y = PC2, group = Treatment), level = 0.8) +
  scale_linetype_manual(values = c(1, 2), guide = FALSE) +
  coord_fixed(1) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
        legend.text.align = 0,
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 14),
        legend.text=element_text(size= 13),
        legend.title = element_text(size = 14))

ggsave('PCA_raw.pdf', width = 13)
ggsave('PCA_raw.jpg', width = 13)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#######################################################################
