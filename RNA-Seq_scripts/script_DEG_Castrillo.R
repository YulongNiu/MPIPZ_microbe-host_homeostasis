
##
## originally by Yulong Niu
## yniu@mpipz.mpg.de
##

##~~~~~~~~~~~~~~~~~~~~~~~~~~~useful funcs~~~~~~~~~~~~~~~~~~~~~~~~
checkZeros <- function(v, threshold) {
    res <- ifelse(sum(v == 0) > threshold, FALSE, TRUE)
    return(res)
}

checkFlg22 <- function(v, threshold) {

    require('magrittr')

    res <- v %>%
    split(rep(1 : 3, each = 6)) %>%
    sapply(checkZeros, threshold) %>%
    all

    return(res)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library('magrittr')
suppressPackageStartupMessages(library('tximport'))
suppressPackageStartupMessages(library('tibble'))
suppressPackageStartupMessages(library('dplyr'))
suppressPackageStartupMessages(library('tidyr'))
suppressPackageStartupMessages(library('readr'))
suppressPackageStartupMessages(library('DESeq2'))
suppressPackageStartupMessages(library('pheatmap'))
suppressPackageStartupMessages(library('goseq'))
suppressPackageStartupMessages(library('foreach'))

anno <- read_csv('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/Ensembl_ath_Anno.csv',
                 col_types = cols(Chromosome = col_character())) %>%
  mutate_all(list(~replace(., is.na(.), '')))

##~~~~~~~~~~~~~~~~~~~~load k alignments~~~~~~~~~~~~~~~~~~~~~~~~~~
wd <- '/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/flg22_crossref/Castrillo_2017/'
setwd(wd)

annoSample <- tibble(condition = paste0(rep(c('Col0', 'flg22'), each = 6),
                                     '_s_',
                                     rep(c('rep1', 'rep2', 'rep3', 'rep1', 'rep2', 'rep3'), each = 2),
                                     rep(c('_batch1', '_batch2'), 6))) %>%
  mutate(sample = paste0(condition, '_ath_kallisto')) %>%
  mutate(batch = rep(1:2, 6)) %>%
  mutate(treatment = rep(c('Col0', 'flg22'), each = 6))

kres <- file.path(wd, annoSample$sample, 'abundance.h5') %>%
  set_names(annoSample$condition) %>%
  tximport(type = 'kallisto', txOut = TRUE)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~normalization~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print('##~~~~~~~~~~~~~~~~~~~~~~~##')
print('step2: normalization counts')
print('##~~~~~~~~~~~~~~~~~~~~~~~##')

condi <- c('Col0', 'flg22')
sampleTable <- data.frame(condition = factor(rep(condi, each = 6), levels = condi))
rownames(sampleTable) <- colnames(kres$counts)

degres <- DESeqDataSetFromTximport(kres, sampleTable, ~condition)

## remove 0|0|x and |0|0|0
degres %<>%
  estimateSizeFactors %>%
  counts(normalized = TRUE) %>%
  apply(1, checkFlg22, 1) %>%
  degres[., ] %>%
  DESeq

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

## ## manual detect surrogate variance
## svnum <- 4
## svseq <- svaseq(dat, mod, mod0, n.sv = svnum)

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

##~~~~~~~~~~~~~~~~~~~~~~~~DEGs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print('##~~~~~~~~~~~~~~~~~~~~~~~##')
print('step3: DEGs')
print('##~~~~~~~~~~~~~~~~~~~~~~~##')

degres$sv1 <- svobj$sv[, 1]
degres$sv2 <- svobj$sv[, 2]
degres$sv3 <- svobj$sv[, 3]
design(degres) <- ~sv1 + sv2 + sv3 + condition

degres <- DESeq(degres)

cond <- list(c('flg22', 'Col0'))

## DEGs
resRaw <- lapply(cond,
                 function(x) {
                   degres %>%
                     results(contrast = c('condition', x)) %T>%
                     summary %>%
                     as_tibble %>%
                     select(pvalue, padj, log2FoldChange) %>%
                     rename_all(.funs =  list(~paste0(paste(x, collapse = '_vs_'), '_', .)))
                 }) %>%
  bind_cols

res <- cbind.data.frame(as.matrix(mcols(degres)[, 1:10]), counts(degres, normalize = TRUE), stringsAsFactors = FALSE) %>%
  rownames_to_column(., var = 'ID') %>%
  as_tibble %>%
  bind_cols(resRaw) %>%
  inner_join(anno, by = 'ID') %>%
  select(ID, Gene : Description, Col0_s_rep1_batch1 : flg22_vs_Col0_log2FoldChange) %>%
  arrange(flg22_vs_Col0_padj)

## selected flg22 DEGs
DEGsflg22 <- res %>%
  filter(abs(flg22_vs_Col0_log2FoldChange) > log2(1.5)) %>%
  filter(flg22_vs_Col0_padj < 0.05)

write_csv(DEGsflg22, 'DEGsflg22_singleend.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~GO enrichment~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print('##~~~~~~~~~~~~~~~~~~~~~~~##')
print('step4: GO enrichment')
print('##~~~~~~~~~~~~~~~~~~~~~~~##')

athGO %<>%
  lapply(function(x) {x[x %in% res$ID]}) %>%
  {
    l <- sapply(., length) > 0
    .[l]
  }

GOMat <- foreach(i = 1:length(athGO), .combine = rbind) %do% {
  eachMat <- cbind(athGO[[i]], names(athGO)[i])
  return(eachMat)
} %>% as.data.frame

## flg22
degVec <- (res$ID %in% DEGsflg22$ID) %>%
  as.integer %>%
  set_names(res$ID)

pdf(NULL)
pwf <- nullp(degVec, bias.data = res$Length)

GOTestWithCat <- goseq(pwf, gene2cat = GOMat, use_genes_without_cat = FALSE) %>%
  as_tibble %>%
  filter(!is.na(ontology))

## print GO
print('Top GO terms flg22 vs Mock')
GOTestWithCat %>%
  head %>%
  as.data.frame %>%
  print

termCat <- c('BP', 'MF', 'CC')
for (i in termCat) {
  write_csv(GOTestWithCat %>% filter(ontology == i),
            paste0(i, 'flg22_singleend.csv') %>% file.path(outdir, .))
}

## MeJA
degVec <- (res$ID %in% DEGsMeJA$ID) %>%
  as.integer %>%
  set_names(res$ID)

pdf(NULL)
pwf <- nullp(degVec, bias.data = res$Length)

GOTestWithCat <- goseq(pwf, gene2cat = GOMat, use_genes_without_cat = FALSE) %>%
  as_tibble %>%
  filter(!is.na(ontology))

## print GO
print('Top GO terms MeJA vs Mock')
GOTestWithCat %>%
  head %>%
  as.data.frame %>%
  print

termCat <- c('BP', 'MF', 'CC')
for (i in termCat) {
  write_csv(GOTestWithCat %>% filter(ontology == i),
            paste0(i, 'MeJA_singleend.csv') %>% file.path(outdir, .))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
