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
wd <- '/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/align_data_1stadd/'
setwd(wd)

annoSample <- read_delim('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/list_samples_1stadd.txt', delim = '\t')

slabel <- annoSample$Anno %>%
  paste0('_ath_kallisto')

files <- file.path(wd, slabel, 'abundance.h5')
names(files) <- annoSample$Anno
idx <- c(1, 11, 21, 31) %>%
  {rep(., 10) + rep(0 : 9, each = 4)}
files %<>% .[idx]
kres <- tximport(files, type = 'kallisto', txOut = TRUE)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~normalization~~~~~~~~~~~~~~~~~~~~~~~
setwd('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/')

## sampleTable
condi <- c('Mock', 'Mock_Flg22', 'HKSynCom33', 'HKSynCom33_Flg22', 'SynCom33', 'SynCom33_Flg22', 'HKSynCom35', 'HKSynCom35_Flg22', 'SynCom35', 'SynCom35_Flg22')
sampleTable <- data.frame(condition = factor(rep(condi, each = 4), levels = condi))
rownames(sampleTable) <- colnames(kres$counts)

degres <- DESeqDataSetFromTximport(kres, sampleTable, ~condition)

## remove 0|0|0|x, 0|0|0|0
degres %<>%
  estimateSizeFactors %>%
  counts(normalized = TRUE) %>%
  apply(1, checkFlg22, 2) %>%
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
ggsave('auto_sv_1stadd.jpg')
ggsave('auto_sv_1stadd.pdf')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~DEG Mock vs. 3 conditions~~~~~~~~~~~~~~~~~
degres$sv1 <- svobj$sv[, 1]
degres$sv2 <- svobj$sv[, 2]
degres$sv3 <- svobj$sv[, 3]
degres$sv4 <- svobj$sv[, 4]
design(degres) <- ~sv1 + sv2 + sv3 + sv4 + condition

cond <- degres %>%
  resultsNames %>%
  str_extract('(?<=condition_).*') %>%
  .[!is.na(.)]

resRaw <- lapply(cond,
                 function(x) {
                   degres %>%
                     results(name = paste0('condition_', x)) %T>%
                     ## lfcShrink(coef = paste0('condition_', x), type = 'apeglm') %T>%
                     summary %>%
                     as_tibble %>%
                     select(pvalue, padj, log2FoldChange) %>%
                     rename_all(.funs = list(~paste0(x, '_', .)))
                 }) %>%
  bind_cols

res <- cbind.data.frame(as.matrix(mcols(degres)[, 1:10]), counts(degres, normalize = TRUE), stringsAsFactors = FALSE) %>%
  rownames_to_column(., var = 'ID') %>%
  as_tibble %>%
  bind_cols(resRaw) %>%
  inner_join(anno, by = 'ID') %>%
  select(ID, Gene : Description, Mock_1 : SynCom35_Flg22_vs_Mock_log2FoldChange) %>%
  arrange(Mock_Flg22_vs_Mock_padj)

write_csv(res, 'eachGroup_vs_Mock_k_1stadd.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

## ## batch correction limma
## rldData %<>% removeBatchEffect(rep(1 : 4, 10) %>% factor)

## ## batch correction sva
## modcombat <- model.matrix(~1, data = sampleTable)
## rldData %<>% ComBat(dat = ., batch = rep(rep(1 : 4, 10) %>% factor) %>% factor, mod = modcombat, par.prior = TRUE, prior.plots = FALSE)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cols <- colData(rld)[, 1] %>% factor(., labels = brewer.pal(10, name = 'Paired'))

## 1 - 2 C
pca <- prcomp(t(rldData))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
percentVar <- round(100 * percentVar)
pca1 <- pca$x[,1]
pca2 <- pca$x[,2]
pcaData <- data.frame(PC1 = pca1, PC2 = pca2, Group = colData(rld)[, 1], ID = rownames(colData(rld)))
ggplot(pcaData, aes(x = PC1, y = PC2, colour = Group, label = ID)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_colour_manual(values = levels(cols)) +
  geom_text_repel(force = 5)
ggsave('PCA_1stadd_sva.pdf', width = 15, height = 12)
ggsave('PCA_1stadd_sva.jpg', width = 15, height = 12)

## 2 - 3 C
pca <- prcomp(t(rldData))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
percentVar <- round(100 * percentVar)
pca1 <- pca$x[,2]
pca2 <- pca$x[,3]
pcaData <- data.frame(PC1 = pca1, PC2 = pca2, Group = colData(rld)[, 1], ID = rownames(colData(rld)))
ggplot(pcaData, aes(x = PC1, y = PC2, colour = Group, label = ID)) +
  geom_point(size = 3) +
  xlab(paste0("PC2: ",percentVar[2],"% variance")) +
  ylab(paste0("PC3: ",percentVar[3],"% variance")) +
  scale_colour_manual(values = levels(cols)) +
  geom_text_repel(force = 5)
ggsave('PCA_1stadd_C23.pdf', width = 15, height = 12)
ggsave('PCA_1stadd_C23.jpg', width = 15, height = 12)


## Mock + Mock_flg22
prefixPCA <- 'M_M22'
sampleIdx <- c(1:8)
colorIdx <- c(1:2)

## 33
prefixPCA <- '33'
sampleIdx <- c(9:24)
colorIdx <- c(3:6)

## 35
prefixPCA <- '35'
sampleIdx <- c(25:40)
colorIdx <- c(7:10)

## flg22
prefixPCA <- 'flg22'
sampleIdx <- c(5:8, 13:16, 21:24, 29:32, 37:40)
colorIdx <- c(2, 4, 6, 8, 10)

## 33 + 33_flg22 + 35 + 35_flg22
prefixPCA <- '33_35'
sampleIdx <- c(17:24, 33:40)
colorIdx <- c(5:6, 9:10)

## mock + mock_flg22 + 33_flg22 + 35_flg22
prefixPCA <- 'like1st'
sampleIdx <- c(1:8, 21:24, 37:40)
colorIdx <- c(1, 2, 6, 10)

rldDataPart <- rldData[, sampleIdx]
pca <- prcomp(t(rldDataPart))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
percentVar <- round(100 * percentVar)
pca1 <- pca$x[, 1]
pca2 <- pca$x[, 2]
pcaData <- data.frame(PC1 = pca1, PC2 = pca2, Group = colData(rld)[sampleIdx, 1], ID = rownames(colData(rld))[sampleIdx])
ggplot(pcaData, aes(x = PC1, y = PC2, colour = Group, label = ID)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_colour_manual(values = levels(cols)[colorIdx]) +
  geom_text_repel(force = 5)
ggsave(paste0('PCA_1stadd_', prefixPCA, '.pdf'), width = 15, height = 12)
ggsave(paste0('PCA_1stadd_', prefixPCA, '.jpg'), width = 15, height = 12)

## 3D plot
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

save(degres, rldData, file = 'degres_condi_Mock_1stadd.RData')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~virus load~~~~~~~~~~~~~~~~~~~~~~~~~~
virusLoad <- read_csv('ath_alignment.csv') %>%
  mutate(K_virus = Kvirus_ath - K_ath)

virusData <- inner_join(virusLoad, pcaData, by = c('sample' = 'ID'))
ggplot(virusData, aes(x = PC1, y = PC2, colour = K_virus, label = sample)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  geom_text_repel(force = 5) +
  scale_colour_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 4, name = 'RdYlBu')))(10)) +
  labs(col = 'Virus load')
ggsave('PCA_1stadd_virusload_sva.pdf', width = 15, height = 12)
ggsave('PCA_1stadd_virusload_sva.jpg', width = 15, height = 12)


hc <- rldData %>%
  cor(method = 'pearson') %>%
  {1 - .} %>%
  as.dist %>%
  hclust(method = 'complete')

hc %>%
  as.dendrogram(method = 'average') %>%
  plot(main = 'Sample Clustering',
       ylab = 'Height')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################
