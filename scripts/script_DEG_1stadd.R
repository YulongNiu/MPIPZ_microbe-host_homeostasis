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
library('tibble')
library('readr')
library('dplyr')
library('stringr')

anno <- read_csv('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/Ensembl_ath_Anno.csv',
                 col_types = cols(Chromosome = col_character())) %>%
  mutate(Gene = Gene %>% {if_else(is.na(.), '', .)}) %>%
  mutate(Description = Description %>% {if_else(is.na(.), '', .)})


##~~~~~~~~~~~~~~~~~~~~load k alignments~~~~~~~~~~~~~~~~~~~~~~~~~~
wd <- '/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/align_data_1stadd'
setwd(wd)

annoSample <- read_delim('/home/Yulong/netscratch/KaWaiFlg22/results/list_samples_1stadd.txt', delim = '\t')

slabel <- annoSample$Anno %>%
  paste0('_ath_kallisto')

files <- file.path(wd, slabel, 'abundance.h5')
names(files) <- annoSample$Anno
idx <- c(1, 11, 21, 31) %>%
  {rep(., 10) + rep(0 : 9, each = 4)}
files %<>% .[idx]
kres <- tximport(files, type = 'kallisto', txOut = TRUE)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~DEG Mock vs. 3 conditions~~~~~~~~~~~~~~~~~
setwd('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/')

## sampleTable
condi <- c('Mock', 'Mock_Flg22', 'HKSynCom33', 'HKSynCom33_Flg22', 'SynCom33', 'SynCom33_Flg22', 'HKSynCom35', 'HKSynCom35_Flg22', 'SynCom35', 'SynCom35_Flg22')
sampleTable <- data.frame(condition = factor(rep(condi, each = 4), levels = condi))
rownames(sampleTable) <- colnames(kres$counts)

degres <- DESeqDataSetFromTximport(kres, sampleTable, ~condition)

## DEGs
degres %<>%
  estimateSizeFactors %>%
  counts(normalized = TRUE) %>%
  apply(1, checkFlg22, 1) %>%
  degres[., ]
## degres <- degres[rowSums(counts(degres)) > 1, ]
save(degres, file = 'degres_condi_Mock_1stadd.RData')
degres <- DESeq(degres)
## resultsNames(degres)

## count transformation
rld <- rlog(degres)
vst <- varianceStabilizingTransformation(degres)
ntd <- normTransform(degres)

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

resRaw <- degres %>%
  results(name = 'condition_Flg22_SynCom35_vs_Flg22_SynCom33') %T>%
  summary %>%
  as_tibble %>%
  select(pvalue, padj, log2FoldChange) %>%
  rename_all(.funs = list(~paste0('Flg22_SynCom35_vs_Flg22_SynCom33_', .)))

res <- cbind.data.frame(as.matrix(mcols(degres)[, 1:10]), assay(ntd), stringsAsFactors = FALSE) %>%
  rownames_to_column(., var = 'ID') %>%
  as_tibble %>%
  bind_cols(resRaw) %>%
  inner_join(anno, by = 'ID') %>%
  select(ID, Gene : Description, Flg22_SynCom33_1 : Flg22_SynCom35_vs_Flg22_SynCom33_log2FoldChange) %>%
  arrange(Flg22_SynCom35_vs_Flg22_SynCom33_padj)

write_csv(res, 'SynCom35_vs_SynCom33_k_full.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PCA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library('directlabels')
library('ggplot2')
library('RColorBrewer')

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## remove low count
thres <- 5
rldData <- assay(rld)

rl <- apply(rldData, 1, function(x){
  return(sum(x > thres) == length(x))
})
rldData %<>% .[rl, ]
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cols <- colData(rld)[, 1] %>% factor(., labels = brewer.pal(10, name = 'Paired'))

pca <- prcomp(t(rldData))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
percentVar <- round(100 * percentVar)
pca1 <- pca$x[,1]
pca2 <- pca$x[,2]
pcaData <- data.frame(PC1 = pca1, PC2 = pca2, Group = colData(rld)[, 1], ID = rownames(colData(rld)))
ggplot(pcaData, aes(x = PC1, y = PC2, colour = Group)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  geom_dl(aes(label = ID, color = Group), method = 'smart.grid') +
  scale_colour_manual(values = levels(cols))
ggsave('PCA_1stadd_thres5.pdf', width = 15, height = 12)
ggsave('PCA_1stadd_thres5.jpg', width = 15, height = 12)


library('rgl')
cols <- colData(rld)[, 1] %>% factor(., labels = brewer.pal(10, name = 'Paired'))
plot3d(pca$x[, 1:3],
       col = cols,
       size = 10,
       xlab = paste0('PC1: ',percentVar[1],'% variance'),
       ylab = paste0('PC2: ',percentVar[2],'% variance'),
       zlab = paste0('PC3: ',percentVar[3],'% variance'))

dir.create('animation_merge')
for (i in 1:360) {
  view3d(userMatrix=rotationMatrix(2*pi * i/360, 0, 1, 0))
  rgl.snapshot(filename=paste('animation_merge/frame-',
                              sprintf('%03d', i), '.jpg', sep=''))}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################
