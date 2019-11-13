###########################DEGs##################################
##~~~~~~~~~~~~~~~~~~~~~~~~~~~useful funcs~~~~~~~~~~~~~~~~~~~~~~~~
checkZeros <- function(v, threshold) {
  res <- ifelse(sum(v == 0) > threshold, FALSE, TRUE)
  return(res)
}

checkFlg22 <- function(v, threshold) {

  require('magrittr')

  res <- v %>%
    split(rep(1 : 4, each = 3)) %>%
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
  mutate_all(list(~replace(., is.na(.), '')))

##~~~~~~~~~~~~~~~~~~~~load k alignments~~~~~~~~~~~~~~~~~~~~~~~~~~
wd <- '/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/align_data'
setwd(wd)

sampleN <- c('Mock', 'Mock_Flg22', 'SynCom33_Flg22', 'SynCom35_Flg22')

slabel <- rep(sampleN, each = 3) %>%
  paste(rep(1:3, 4), sep = '_') %>%
  paste0('_ath_kallisto')

files <- file.path(wd, slabel, 'abundance.h5')
names(files) <- rep(sampleN, each = 3) %>%
  paste(rep(1:3, 4), sep = '_')
kres <- tximport(files, type = 'kallisto', txOut = TRUE)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~normalization~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/')

condi <- c('Mock', 'Mock_Flg22', 'SynCom33_Flg22', 'SynCom35_Flg22')

## sampleTable
sampleTable <- data.frame(condition = factor(rep(condi, each = 3))) %>%
  set_rownames(colnames(kres$counts))

sampleTable$condition %<>% relevel(ref = 'Mock')

degres <- DESeqDataSetFromTximport(kres, sampleTable, ~condition)

## remove 0|0|x, 0|0|0
degres %<>%
  estimateSizeFactors %>%
  counts(normalized = TRUE) %>%
  apply(1, checkFlg22, 1) %>%
  degres[., ]

## save(degres, file = 'degres_condi_Mock.RData')
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
ggsave('auto_sv.jpg')
ggsave('auto_sv.pdf')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~DEG Mock vs. 3 conditions~~~~~~~~~~~~~~~~~
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
  select(ID, Gene : Description, Mock_1 : Flg22_SynCom35_vs_Mock_log2FoldChange) %>%
  arrange(Flg22_vs_Mock_padj)

write_csv(res, 'eachGroup_vs_Mock_k.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PCA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library('directlabels')
library('ggplot2')

pca <- prcomp(t(assay(rld)))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
percentVar <- round(100 * percentVar)
pca1 <- pca$x[,1]
pca2 <- pca$x[,2]
pcaData <- data.frame(PC1 = pca1, PC2 = pca2, Group = colData(rld)[, 1], ID = rownames(colData(rld)))
cairo_pdf('PCA.pdf', width = 12)
ggplot(pcaData, aes(x = PC1, y = PC2, colour = Group)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  geom_dl(aes(label = ID, color = Group), method = 'smart.grid')
dev.off()
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################
