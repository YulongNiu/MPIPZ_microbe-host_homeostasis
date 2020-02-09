###########################DEGs##################################
##~~~~~~~~~~~~~~~~~~~~~~~~~~~useful funcs~~~~~~~~~~~~~~~~~~~~~~~~
checkZeros <- function(v, threshold) {
  res <- ifelse(sum(v == 0) > threshold, FALSE, TRUE)
  return(res)
}

checkFlg22 <- function(v, threshold) {

  require('magrittr')

  res <- v %>%
    split(rep(1 : 3, each = 3)) %>%
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

anno <- read_csv('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/Ensembl_ath_Anno.csv',
                 col_types = cols(Chromosome = col_character())) %>%
  mutate(Gene = Gene %>% {if_else(is.na(.), '', .)}) %>%
  mutate(Description = Description %>% {if_else(is.na(.), '', .)})


##~~~~~~~~~~~~~~~~~~~~load k alignments~~~~~~~~~~~~~~~~~~~~~~~~~~
wd <- '/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/align_data_soil'
setwd(wd)

annoSample <- read_csv('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/Ka-Wai_soil.csv')

slabel <- annoSample$Anno %>%
  paste0('_ath_kallisto')

files <- file.path(wd, slabel, 'abundance.h5')
names(files) <- annoSample$Anno
idx <- c(1, 4, 7, 2, 5, 8, 3, 6, 9)
files %<>% .[idx]
kres <- tximport(files, type = 'kallisto', txOut = TRUE)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~DEG Mock vs. 3 conditions~~~~~~~~~~~~~~~~~
setwd('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/')

## sampleTable
condi <- c('Mock', 'SynCom33', 'SynCom35')
sampleTable <- data.frame(condition = factor(rep(condi, each = 3), levels = condi))
rownames(sampleTable) <- colnames(kres$counts)

sampleTable$condition %<>% relevel(ref = 'Mock')

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

##~~~~~~~~~~~~~~~~~~~~~~~~~~hidden batch effect~~~~~~~~~~~~~~~~~~~~~
library('sva')
library('ggplot2')

dat <- rld %>%
  assay %>%
  {.[rowMeans(.) > 1, ]}
mod <- model.matrix(~ condition, colData(degres))
mod0 <- model.matrix(~ 1, colData(degres))

## manual detect surrogate variance
svnum <- 4
svseq <- svaseq(dat, mod, mod0, n.sv = svnum)

## auto detect sv
svobj <- sva(dat, mod, mod0)
svnum <- svobj$sv %>% ncol

svobj$sv %>%
  set_colnames(paste0('sv', seq_len(svnum))) %>%
  as_tibble %>%
  gather(key = 'sv', value = 'value') %>%
  mutate(condition = colData(degres) %>%
           .$condition %>%
           rep(svnum) %>%
           as.character,
         sample = rep(colnames(degres), svnum)) %>%
  mutate(group = paste(sv, condition, sep = '_')) %>%
  ggplot(aes(sample, value, colour = sv, group = group)) +
  geom_point() +
  geom_line() +
  theme(axis.text.x = element_text(angle = 90))
ggsave('auto_sv_soil.jpg')
ggsave('auto_sv_soil.pdf')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cond <- degres %>%
  resultsNames %>%
  str_extract('(?<=condition_).*') %>%
  .[!is.na(.)]

resRaw <- lapply(cond,
                 function(x) {
                   degres %>%
                     results(name = paste0('condition_', x)) %T>%
                     summary %>%
                     as_tibble %>%
                     select(pvalue, padj, log2FoldChange) %>%
                     rename_all(.funs = list(~paste0(x, '_', .)))
                 }) %>%
  bind_cols


res <- cbind.data.frame(as.matrix(mcols(degres)[, 1:10]), assay(ntd), stringsAsFactors = FALSE) %>%
  rownames_to_column(., var = 'ID') %>%
  as_tibble %>%
  bind_cols(resRaw) %>%
  inner_join(anno, by = 'ID') %>%
  select(ID, Gene : Description, Mock_1 : SynCom35_Flg22_vs_Mock_log2FoldChange) %>%
  arrange(Mock_Flg22_vs_Mock_padj)

write_csv(res, 'eachGroup_vs_Mock_k_1stadd.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


##~~~~~~~~~~~~~~~~~~~~~~~~DEG flg22_SynCom vs flg22~~~~~~~~~~~~~~~~~
setwd('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/')

## sampleTable
sampleTable <- data.frame(condition = factor(rep(c('Mock', 'Flg22', 'Flg22_SynCom33', 'Flg22_SynCom35'), each = 3)), process = factor(rep(c('N', 'Y', 'Y', 'Y'), each = 3)))
sampleTable$condition %<>% relevel(ref = 'Flg22')
rownames(sampleTable) <- colnames(kres$counts)

degres <- DESeqDataSetFromTximport(kres,
                                   sampleTable,
                                   ~condition)

## DEGs role out two zeros in one group
## degres <- degres[rowSums(counts(degres)) > 1, ]
degres %<>%
  estimateSizeFactors %>%
  counts(normalized = TRUE) %>%
  apply(1, checkFlg22, 1) %>%
  degres[., ]
degres <- DESeq(degres)

## count transformation
rld <- rlog(degres)
vst <- varianceStabilizingTransformation(degres)
ntd <- normTransform(degres)

cond <- degres %>%
  resultsNames %>%
  str_extract('(?<=condition_).*') %>%
  .[!is.na(.)] %>%
  `[`(1:2)

resRaw <- lapply(cond,
                 function(x) {
                   degres %>%
                     results(name = paste0('condition_', x)) %T>%
                     summary %>%
                     as_tibble %>%
                     select(pvalue, padj, log2FoldChange) %>%
                     rename_all(.funs = list(~paste0(x, '_', .)))
                 }) %>%
  bind_cols

res <- cbind.data.frame(as.matrix(mcols(degres)[, 1:10]), assay(ntd), stringsAsFactors = FALSE) %>%
  rownames_to_column(., var = 'ID') %>%
  as_tibble %>%
  bind_cols(resRaw) %>%
  inner_join(anno, by = 'ID') %>%
  select(ID, Gene : Description, Flg22_1 : Flg22_SynCom35_vs_Flg22_log2FoldChange) %>%
  arrange(Flg22_SynCom33_vs_Flg22_padj)

write_csv(res, 'SynCom_vs_flg22_k.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~SynCom35 vs. SynCom35~~~~~~~~~~~~~~~~~~
setwd('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/')

## sampleTable
sampleTable <- data.frame(condition = factor(rep(c('Mock', 'Flg22', 'Flg22_SynCom33', 'Flg22_SynCom35'), each = 3)), process = factor(rep(c('N', 'Y', 'Y', 'Y'), each = 3)))
sampleTable$condition %<>% relevel(ref = 'Flg22_SynCom33')
rownames(sampleTable) <- colnames(kres$counts)

degres <- DESeqDataSetFromTximport(kres,
                                   sampleTable,
                                   ~condition)

## DEGs role out two zeros in one group
## degres <- degres[rowSums(counts(degres)) > 1, ]
degres %<>%
  estimateSizeFactors %>%
  counts(normalized = TRUE) %>%
  apply(1, checkFlg22, 1) %>%
  degres[., ]
degres <- DESeq(degres)

## count transformation
rld <- rlog(degres)
vst <- varianceStabilizingTransformation(degres)
ntd <- normTransform(degres)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~DEG Mock vs. 3 conditions~~~~~~~~~~~~~~~~~
degres$sv1 <- svobj$sv[, 1]
degres$sv2 <- svobj$sv[, 2]
design(degres) <- ~sv1 + sv2 + condition

degres <- DESeq(degres)

cond <- list(c('SynCom33', 'Mock'),
             c('SynCom35', 'Mock'),
             c('SynCom33', 'SynCom35'))

resRaw <- lapply(cond,
                 function(x) {
                   degres %>%
                     results(contrast = c('condition', x)) %T>%
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
  select(ID, Gene : Description, Mock_1 : SynCom33_vs_SynCom35_log2FoldChange) %>%
  arrange(SynCom33_vs_Mock_padj)

write_csv(res, 'eachGroup_vs_Mock_k_1stadd.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PCA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library('ggrepel')
library('ggplot2')
library('RColorBrewer')
library('limma')
library('sva')

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## remove low count
dat <- rld %>%
  assay %>%
  {.[rowMeans(.) > 1, ]}

group <- sampleTable$condition
design <- model.matrix(~ group)
rldData <- dat %>%
  removeBatchEffect(covariates = svobj$sv,
                    design = design)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pca <- prcomp(t(rldData))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
percentVar <- round(100 * percentVar)
pca1 <- pca$x[,1]
pca2 <- pca$x[,2]
pcaData <- data.frame(PC1 = pca1, PC2 = pca2, Group = colData(rld)[, 1], ID = rownames(colData(rld))) %>%
  mutate(SynCom = rep(c('Mock', 'Non-suppressive', 'Suppressive'), c(3, 3, 3))) %>%
  mutate(Treatment = rep(c('Mock', 'SynCom'), c(3, 6)) %>% factor) %>%
  mutate(flg22 = c('without') %>% rep(9) %>% factor)

ggplot(pcaData, aes(x = PC1, y = PC2, colour = SynCom)) +
  geom_point(aes(shape = Treatment), size = 4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_colour_manual(values = c('#000000', '#377eb8', '#e41a1c')) +
  scale_shape_manual(values = c(1, 15),
                     name = 'Experimental\nConditions') +
  stat_ellipse(aes(x = PC1, y = PC2, group = Group, linetype = flg22), type = 't', level = 0.7) +
  scale_linetype_manual(values = c(1, 2), guide = FALSE) +
  coord_fixed(1) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
        legend.text.align = 0,
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 14),
        legend.text=element_text(size= 13),
        legend.title = element_text(size = 14))

ggsave('PCA_soil_sva.pdf', width = 13)
ggsave('PCA_soil_sva.jpg', width = 13)

save(degres, rldData, file = 'degres_condi_Mock.RData')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################
