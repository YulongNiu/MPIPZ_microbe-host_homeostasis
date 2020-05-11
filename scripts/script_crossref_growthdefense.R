
## originally by Yulong Niu
## yulong.niu@hotmail.com

##########################cross ref growth defense######################
library('tidyverse')

setwd('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/removeZero/growth_defense/')

## 1 AT4G34390.1 AT4G34390 XLG2
## 2 AT2G17720.1 AT2G17720 P4H5
## 3 AT3G53180.1 AT3G53180 AT3G53180
## 4 AT5G51060.1 AT5G51060 RBOHC
gdGenes <- c('AT4G34390', 'AT2G17720', 'AT3G53180', 'AT5G51060')

##~~~~~~~~~~~~~~~~~~~~~~~~~~Col0 agar plate~~~~~~~~~~~~~~~~~~~~~~~~~~~
load('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/removeZero/degres_condi_Mock_1stadd.RData')

rawC <- rldData[str_detect(rownames(rldData), paste(gdGenes, collapse = '|')), ] %>%
  as.data.frame %>%
  rownames_to_column('ID') %>%
  as_tibble %>%
  dplyr::select(matches('Mock_\\d|HKSynCom33_\\d|HKSynCom35_\\d'), matches('Mock_Flg22_\\d|HKSynCom33_Flg22_\\d|HKSynCom35_Flg22_\\d'), everything())

Col0Agar <- rawC %>%
  select(contains('_')) %>%
  t %>%
  scale %>%
  t %>%
  as_tibble %>%
  bind_cols(rawC %>% select(ID))

Col0Agar %>%
  gather(-ID, key = 'Sample', value = 'ScaleCounts') %>%
  mutate(Conditions = case_when(
           str_detect(Sample, '^Mock_\\d') ~ 'Mock',
           str_detect(Sample, '^Mock_Flg22') ~ 'Mock_flg22',
           str_detect(Sample, '^HKSynCom33_\\d') ~ 'HK_nonsupp',
           str_detect(Sample, '^HKSynCom33_Flg22') ~ 'HK_nonsupp_flg22',
           str_detect(Sample, '^HKSynCom35_\\d') ~ 'HK_supp',
           str_detect(Sample, '^HKSynCom35_Flg22') ~ 'HK_supp_flg22',
           str_detect(Sample, '^SynCom33_\\d') ~ 'Nonsupp',
           str_detect(Sample, '^SynCom33_Flg22') ~ 'Nonsupp_flg22',
           str_detect(Sample, '^SynCom35_\\d') ~ 'Supp',
           str_detect(Sample, '^SynCom35_Flg22') ~ 'Supp_flg22'
         )) %>%
  mutate(Group = case_when(
           str_detect(Conditions, '^Mock|^HK') ~ 'Mock',
           str_detect(Conditions, '^Nonsupp') ~ 'Nonsupp',
           str_detect(Conditions, '^Supp') ~ 'Supp'
         )) %>%
  mutate(Conditions = Conditions %>% factor(levels = c('Mock', 'HK_nonsupp',  'HK_supp', 'Mock_flg22', 'HK_nonsupp_flg22', 'HK_supp_flg22', 'Nonsupp', 'Nonsupp_flg22', 'Supp', 'Supp_flg22'))) %>%
  ggplot(aes(x = Group, y = ScaleCounts, fill = Conditions)) +
  geom_boxplot(position = position_dodge2(preserve = 'single')) +
  ## geom_boxplot(position = position_dodge(0.8)) +
  ## geom_dotplot(binaxis = 'y', stackdir = 'center', position = position_dodge(0.8), dotsize = 0.5) +
  scale_fill_manual(values = c(rep(NA, 6), rep('#377eb8', 2), rep('#e41a1c', 2))) +
  facet_wrap(~ ID, ncol = 2) +
  ylab('Scaled counts') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90),
        plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
        legend.text.align = 0,
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 14),
        legend.text=element_text(size= 13),
        legend.title = element_text(size = 14))

ggsave(file = 'Ka-Wai_Col0_Agar.pdf')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~pWER agar~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/removeZero/degres_condi_Mock.RData')

rawC <- rldData[str_detect(rownames(rldData), paste(gdGenes, collapse = '|')), ] %>%
  as.data.frame %>%
  rownames_to_column('ID') %>%
  as_tibble

pWERAgar <- rawC %>%
  select(contains('_')) %>%
  t %>%
  scale %>%
  t %>%
  as_tibble %>%
  bind_cols(rawC %>% select(ID))

pWERAgar %>%
  gather(-ID, key = 'Sample', value = 'ScaleCounts') %>%
  mutate(Conditions = case_when(
           str_detect(Sample, '^Mock_\\d') ~ 'Mock',
           str_detect(Sample, '^Mock_Flg22') ~ 'Mock_flg22',
           str_detect(Sample, '^SynCom33_Flg22') ~ 'Nonsupp_flg22',
           str_detect(Sample, '^SynCom35_Flg22') ~ 'Supp_flg22'
         )) %>%
  mutate(Group = case_when(
           str_detect(Conditions, '^Mock|^HK') ~ 'Mock',
           str_detect(Conditions, '^Nonsupp') ~ 'Nonsupp',
           str_detect(Conditions, '^Supp') ~ 'Supp'
         )) %>%
  mutate(Conditions = Conditions %>% factor(levels = c('Mock', 'Mock_flg22', 'Nonsupp_flg22', 'Supp_flg22'))) %>%
  ggplot(aes(x = Group, y = ScaleCounts, fill = Conditions)) +
  geom_boxplot(position = position_dodge2(preserve = 'single')) +
  ## geom_boxplot(position = position_dodge(0.8)) +
  ## geom_dotplot(binaxis = 'y', stackdir = 'center', position = position_dodge(0.8), dotsize = 0.5) +
  scale_fill_manual(values = c(rep(NA, 2), rep('#377eb8', 1), rep('#e41a1c', 1))) +
  facet_wrap(~ ID, ncol = 2) +
  ylab('Scaled counts') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90),
        plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
        legend.text.align = 0,
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 14),
        legend.text=element_text(size= 13),
        legend.title = element_text(size = 14))

ggsave(file = 'Ka-Wai_pWER_Agar.pdf')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~CJ iron ion starvation~~~~~~~~~~~~~~~~~~
load('/extDisk1/RESEARCH/MPIPZ_CJ_RNASeq/results/eachGroup_mergeDay8.RData')

rawC <- rldData[str_detect(rownames(rldData), paste(gdGenes, collapse = '|')), ] %>%
  as.data.frame %>%
  rownames_to_column('ID') %>%
  as_tibble

iron <- rawC %>%
  select(contains('_')) %>%
  t %>%
  scale %>%
  t %>%
  as_tibble %>%
  bind_cols(rawC %>% select(ID))

iron %>%
  gather(-ID, key = 'Sample', value = 'ScaleCounts') %>%
  mutate(Conditions = case_when(
           str_detect(Sample, 'Col0_FeCl3_HK') ~ 'Col0_FeCl3_HK',
           str_detect(Sample, 'Col0_FeCl3_Live') ~ 'Col0_FeCl3_Live',
           str_detect(Sample, 'Col0_FeEDTA_HK') ~ 'Col0_FeEDTA_HK',
           str_detect(Sample, 'Col0_FeEDTA_Live') ~ 'Col0_FeEDTA_Live',
           str_detect(Sample, 'f6h1_FeCl3_HK') ~ 'f6h1_FeCl3_HK',
           str_detect(Sample, 'f6h1_FeCl3_Live') ~ 'f6h1_FeCl3_Live',
           str_detect(Sample, 'f6h1_FeEDTA_HK') ~ 'f6h1_FeEDTA_HK',
           str_detect(Sample, 'f6h1_FeEDTA_Live') ~ 'f6h1_FeEDTA_Live'
         )) %>%
  mutate(Group = case_when(
           str_detect(Conditions, '^Col0') ~ 'Col0',
           str_detect(Conditions, '^f6h1') ~ 'f6h1'
         )) %>%
  mutate(Conditions = Conditions %>% factor(levels = c('Col0_FeCl3_HK', 'Col0_FeCl3_Live', 'Col0_FeEDTA_HK', 'Col0_FeEDTA_Live', 'f6h1_FeCl3_HK', 'f6h1_FeCl3_Live', 'f6h1_FeEDTA_HK', 'f6h1_FeEDTA_Live'))) %>%
  ggplot(aes(x = Group, y = ScaleCounts, fill = Conditions)) +
  geom_boxplot(position = position_dodge2(preserve = 'single')) +
  ## geom_boxplot(position = position_dodge(0.8)) +
  ## geom_dotplot(binaxis = 'y', stackdir = 'center', position = position_dodge(0.8), dotsize = 0.5) +
  scale_fill_manual(values = c(NA, '#e41a1c', NA, '#377eb8') %>% rep(2)) +
  facet_wrap(~ ID, ncol = 2) +
  ylab('Scaled counts') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90),
        plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
        legend.text.align = 0,
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 14),
        legend.text=element_text(size= 13),
        legend.title = element_text(size = 14))

ggsave(file = 'CJ_iron.pdf')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Kathrin host microbia~~~~~~~~~~~~~~~~~~
load('/extDisk1/RESEARCH/MPIPZ_Kathrin_Persistence_RNASeq/results/degres_condi_Mock_ath.RData')

rawC <- rldData[str_detect(rownames(rldData), paste(gdGenes, collapse = '|')), ] %>%
  as.data.frame %>%
  rownames_to_column('ID') %>%
  as_tibble

hostmicro <- rawC %>%
  select(contains('_')) %>%
  t %>%
  scale %>%
  t %>%
  as_tibble %>%
  bind_cols(rawC %>% select(ID))

hostmicro %>%
  gather(-ID, key = 'Sample', value = 'ScaleCounts') %>%
  mutate(Conditions = case_when(
           str_detect(Sample, '^C_fSC') ~ 'fullSC',
           str_detect(Sample, '^C_AtSC') ~ 'AtSC',
           str_detect(Sample, '^C_LjSC') ~ 'LjSC',
           str_detect(Sample, '^C_mock') ~ 'Mock'
         )) %>%
  mutate(Group = case_when(
           str_detect(Conditions, '^Mock') ~ 'Mock',
           str_detect(Conditions, 'fullSC|AtSC|LjSC') ~ 'SynCom'
         )) %>%
  mutate(Conditions = Conditions %>% factor(levels = c('Mock', 'AtSC', 'LjSC', 'fullSC'))) %>%
  ggplot(aes(x = Group, y = ScaleCounts, fill = Conditions)) +
  geom_boxplot(position = position_dodge2(preserve = 'single')) +
  ## geom_boxplot(position = position_dodge(0.8)) +
  ## geom_dotplot(binaxis = 'y', stackdir = 'center', position = position_dodge(0.8), dotsize = 0.5) +
  scale_fill_manual(values = c(NA, '#377eb8', '#e41a1c', '#984ea3')) +
  facet_wrap(~ ID, ncol = 2) +
  ylab('Scaled counts') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90),
        plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
        legend.text.align = 0,
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 14),
        legend.text=element_text(size= 13),
        legend.title = element_text(size = 14))

ggsave(file = 'Kathrin_hostmicrobia.pdf')
########################################################################
