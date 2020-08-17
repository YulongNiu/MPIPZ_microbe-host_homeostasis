### Scripts originally created by Ruben Garrido-Oter and modified by Ka-Wai Ma
# The scipt can be used to generate the PCA and CPcoA plots. Change the name of
# the files and samples e.g. strain ID accordingly for all 3 experiments

options(warn=-1)

# cleanup

rm(list=ls())

# load plotting functions
source("plotting_functions.R")
source("cpcoa.func.R")

if (!require(tidyverse)) install.packages("tidyverse")
if (!require(FSA)) install.packages("FSA")
if (!require(multcomp)) install.packages("multcomp")

library(tidyverse)
library(multcomp)
library(FSA)

# directories

results.dir <- "../results/"
data.dir <- "../data/"
figures.dir <- "../figures/"

# files

design.file <- paste(data.dir, "pwer_001_design.txt", sep="")
otu_table.file <- paste(results.dir, "otu_table_pwer1.txt", sep="")
otu_table_unfiltered.file <- paste(results.dir, "otu_table_unfiltered_pwer1.txt", sep="")
taxonomy.file <- paste(data.dir, "pwer_001_strains.txt", sep="")

# load data

design <- read.table(design.file, header=T, sep="\t")
otu_table <- read.table(otu_table.file, sep="\t", header=T, check.names=F)
otu_table_unfiltered <- read.table(otu_table_unfiltered.file, sep="\t", header=T, check.names=F)
taxonomy <- read.table(taxonomy.file, sep="\t", header=T, check.names=F)

# re-order data matrices

idx <- design$SampleID %in% colnames(otu_table)
design <- design[idx, ]
idx <- match(design$SampleID, colnames(otu_table))
otu_table <- otu_table[, idx]
idx <- match(design$SampleID, colnames(otu_table_unfiltered))
otu_table_unfiltered <- otu_table_unfiltered[, idx]

# normalize otu tables

design$depth <- colSums(otu_table)
design$depth_unfiltered <- colSums(otu_table_unfiltered)
otu_table_norm <- apply(otu_table, 2, function(x) x/sum(x))
otu_table_unfiltered_norm <- apply(otu_table_unfiltered, 2, function(x) x/sum(x))

idx <- rownames(otu_table_unfiltered_norm) %in% c("good", "noisy", "chimera", "other")
design <- cbind(design, t(otu_table_unfiltered_norm[idx, ]))
design$exact <- 1-rowSums(t(otu_table_unfiltered_norm[!idx, ]))

idx <- rownames(otu_table_norm) == "ecoli_dh5a"
ecoli <- otu_table_norm[idx, ]
otu_table <- otu_table[!idx,  ]
otu_table_norm <- apply(otu_table, 2, function(x) x/sum(x))

design <- cbind(design, ecoli)

design$ecoli_density <- design$ecoli / design$root.weight


idx <- rownames(otu_table_norm) %in% taxonomy$strain[taxonomy$syncom=="NS"]
aggRANS <- colSums(otu_table_norm[idx, ])
design <- cbind(design, aggRANS)
idx <- rownames(otu_table_norm) %in% taxonomy$strain[taxonomy$syncom=="S"]
aggRAS <- colSums(otu_table_norm[idx, ])
design <- cbind(design, aggRAS)


design$cond <- paste(design$SynCom, design$flg22, design$genotype)

### compute beta diversity
colors <- data.frame(group=c("input", "agar", "agar_plant","root"),
                     color=c(c_black, c_grey, c_dark_brown, c_orange))

# PCoA Bray-Curtis
# subset samples of interest
# compartment-specific effect of samples inoculated with different SynComs
# we will use all data point except for replicate c in experiment 1 due to
# potential contamination or PCR error.

# SynCom NS
idx <- design$compartment %in% c("input", "root", "agar", "agar_plant") &
  design$SynCom %in% c("NS") &
  !design$biological.replicate=="c" &  
  TRUE

design_fullNS <- design[idx, ]
otu_table_norm_fullNS <- otu_table_norm[, idx]
otu_table_unfiltered_norm_NS_full <- otu_table_unfiltered_norm[, idx]

bray_curtis <- vegdist(t(otu_table_norm_fullNS), method="bray")

k <- 2
pcoa <- cmdscale(bray_curtis, k=k, eig=T)
points <- pcoa$points
eig <- pcoa$eig
points <- as.data.frame(points)
colnames(points) <- c("x", "y")

points_fullNS <- cbind(points, design[match(rownames(points), design$SampleID), ])

colors_pcoa <- colors[colors$group %in% points_fullNS$compartment, ]
points_fullNS$compartment <- factor(points_fullNS$compartment, levels=colors_pcoa$group)

# SynCom S
idx <- design$compartment %in% c("input", "root", "agar", "agar_plant") &
  design$SynCom %in% c("S") &
  !design$biological.replicate =="c" &  
    TRUE

design_fullS <- design[idx, ]
otu_table_norm_fullS <- otu_table_norm[, idx]
otu_table_unfiltered_norm_S_full <- otu_table_unfiltered_norm[, idx]

bray_curtis <- vegdist(t(otu_table_norm_fullS), method="bray")
k <- 2
pcoa <- cmdscale(bray_curtis, k=k, eig=T)
points <- pcoa$points
eig <- pcoa$eig
points <- as.data.frame(points)
colnames(points) <- c("x", "y")

points_fullS <- cbind(points, design[match(rownames(points), design$SampleID), ])

colors_pcoa <- colors[colors$group %in% points_fullS$compartment, ]
points_fullS$compartment <- factor(points_fullS$compartment, levels=colors_pcoa$group)

# SynCom NS+S
idx <- design$compartment %in% c("input", "root", "agar", "agar_plant") &
  design$SynCom %in% c("NS+S") &
  !design$biological.replicate=="c" &  
    TRUE

design_fullNSS <- design[idx, ]
otu_table_norm_fullNSS <- otu_table_norm[, idx]
otu_table_unfiltered_norm_NSS_full <- otu_table_unfiltered_norm[, idx]

bray_curtis <- vegdist(t(otu_table_norm_fullNSS), method="bray")
k <- 2
pcoa <- cmdscale(bray_curtis, k=k, eig=T)
points <- pcoa$points
eig <- pcoa$eig
points <- as.data.frame(points)
colnames(points) <- c("x", "y")

points_fullNSS <- cbind(points, design[match(rownames(points), design$SampleID), ])

colors_pcoa <- colors[colors$group %in% points_fullNSS$compartment, ]
points_fullNSS$compartment <- factor(points_fullNSS$compartment, levels=colors_pcoa$group)

### plot PCo 1 and 2

ggplot(points_fullNS, aes(x=x, y=y, color=compartment)) +
  geom_point(alpha=.7, size=4.0) +
  scale_colour_manual(values=as.character(colors_pcoa$color)) +
  theme(axis.text.x=element_text(size=12)) +
  theme(axis.text.y=element_text(size=12)) +
  theme(axis.title=element_text(size=13)) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
  main_theme +
  theme(legend.position="top") +
  ggsave(paste(figures.dir, "PCoA_NS.pdf", sep=""), width=7, height=5)
#  ggsave(paste(figures.dir, "PCoA_NS.tiff", sep=""), width=7, height=5)

ggplot(points_fullS, aes(x=x, y=y, color=compartment)) +
  geom_point(alpha=.7, size=4.0) +
  scale_colour_manual(values=as.character(colors_pcoa$color)) +
  theme(axis.text.x=element_text(size=12)) +
  theme(axis.text.y=element_text(size=12)) +
  theme(axis.title=element_text(size=13)) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
  main_theme +
  theme(legend.position="top") +
  ggsave(paste(figures.dir, "PCoA_S.pdf", sep=""), width=7, height=5)
#  ggsave(paste(figures.dir, "PCoA_S.tiff", sep=""), width=7, height=5)

ggplot(points_fullNSS, aes(x=x, y=y, color=compartment)) +
  geom_point(alpha=.7, size=4.0) +
  scale_colour_manual(values=as.character(colors_pcoa$color)) +
  theme(axis.text.x=element_text(size=12)) +
  theme(axis.text.y=element_text(size=12)) +
  theme(axis.title=element_text(size=13)) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
  main_theme +
  theme(legend.position="top") +
  ggsave(paste(figures.dir, "PCoA_NSS.pdf", sep=""), width=7, height=5)
#  ggsave(paste(figures.dir, "PCoA_NSS.tiff", sep=""), width=7, height=5)


### CPCoA Bray-Curtis of samples conditioned by technical factors
colors <- data.frame(group=c(
  "NS FALSE pWER::FLS2-GFP",
  "NS TRUE pWER::FLS2-GFP",
  "NS FALSE Col-0",
  "NS TRUE Col-0",
  "S FALSE pWER::FLS2-GFP",
  "S TRUE pWER::FLS2-GFP",
  "S FALSE Col-0",
  "S TRUE Col-0",
  "NS+S FALSE pWER::FLS2-GFP",
  "NS+S TRUE pWER::FLS2-GFP",
  "NS+S FALSE Col-0",
  "NS+S TRUE Col-0"
),
color=c(
  "#f8766d",
  "#f8766d",
  "#f8766d",
  "#f8766d",
  "#00bfc4",
  "#00bfc4",
  "#00bfc4",
  "#00bfc4",
  "#800080",
  "#800080",
  "#800080",
  "#800080"
))

shapes <- data.frame(group=c(
  "NS FALSE pWER::FLS2-GFP",
  "NS TRUE pWER::FLS2-GFP",
  "NS FALSE Col-0",
  "NS TRUE Col-0",
  "S FALSE pWER::FLS2-GFP",
  "S TRUE pWER::FLS2-GFP",
  "S FALSE Col-0",
  "S TRUE Col-0",
  "NS+S FALSE pWER::FLS2-GFP",
  "NS+S TRUE pWER::FLS2-GFP",
  "NS+S FALSE Col-0",
  "NS+S TRUE Col-0"
),
shape=c(
  2,
  17,
  1,
  16,
  2,
  17,
  1,
  16,
  2,
  17,
  1,
  16
))

sqrt_transform <- T

# CPcoA for root samples in different genotypes 
idx <- design$compartment %in% c("root") &  
  design$genotype %in% c("pWER::FLS2-GFP") &
  #design$genotype %in% c("Col-0") &
  !design$biological.replicate %in% c("c") &  
  T


d <- design[idx, ]
idx <- colnames(otu_table_norm) %in% d$SampleID

bray_curtis <- vegdist(t(otu_table_norm[, idx]), method="bray")

capscale.gen <- capscale(bray_curtis ~ flg22 * SynCom + Condition(biological.replicate * technical.replicate), data=d, add=F, sqrt.dist=sqrt_transform)

perm_anova.gen <- anova.cca(capscale.gen)

# generate variability tables and calculate confidence intervals for the variance
var_tbl.gen <- variability_table(capscale.gen)

eig <- capscale.gen$CCA$eig

variance <- capscale.gen$CCA$tot.chi / capscale.gen$tot.chi

variance <- var_tbl.gen["constrained", "proportion"]
p.val <- perm_anova.gen[1, 4]

points_cpcoa <- capscale.gen$CCA$wa[, 1:2]
colnames(points_cpcoa) <- c("x", "y")
points_cpcoa <- cbind(points_cpcoa, d[match(rownames(points_cpcoa), d$SampleID), ])

# plot CPCo 1 and 2

colors_cpcoa <- colors[colors$group %in% points_cpcoa$cond, ]
points_cpcoa$cond1 <- factor(points_cpcoa$cond, levels=colors_cpcoa$group)

shapes_cpcoa <- shapes[shapes$group %in% points_cpcoa$cond, ]
points_cpcoa$cond2 <- factor(points_cpcoa$cond, levels=shapes_cpcoa$group)


ggplot(points_cpcoa, aes(x=x, y=y, color=cond1, shape=cond2)) +
  geom_point(alpha=.7, size=4.0) +
  stat_ellipse(aes(color=cond1), linetype=2, type="t", size=1.0, alpha=.5) +
  scale_colour_manual(values=as.character(colors_cpcoa$color)) +
  scale_shape_manual(values=shapes_cpcoa$shape) +
  theme(axis.text.x=element_text(size=12)) +
  theme(axis.text.y=element_text(size=12)) +
  theme(axis.title=element_text(size=13)) +
  labs(x=paste("CPCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("CPCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) +
  ggtitle("CPCoA of root samples in pwer", paste(format(100 * variance, digits=3), " % of variance; P=", format(p.val, digits=2), sep="")) +
  main_theme +
  theme(legend.position="bottom")+
  ggsave(paste(figures.dir, "CPCoA_root_pwert.pdf", sep=""),  height=6, width=5)
#  ggsave(paste(figures.dir, "CPCoA_root_pwert.tiff", sep=""),  height=6, width=5)

# Repeat the above cPCoA for samples in different geneotype including Col-0

# Permutation analysis for samples inoculated by different SynCom e.g. NS, S and NS+S

idx <- design$compartment %in% c("root") &  
  design$genotype %in% c("pWER::FLS2-GFP") &
  design$SynCom %in% c("NS+S") &
  !design$biological.replicate %in% c("c") &  
  T

d <- design[idx, ]
idx <- colnames(otu_table_norm) %in% d$SampleID

bray_curtis <- vegdist(t(otu_table_norm[, idx]), method="bray")

capscale.gen <- capscale(bray_curtis ~ flg22 * SynCom + Condition(biological.replicate * technical.replicate), data=d, add=F, sqrt.dist=sqrt_transform)

perm_anova.gen <- anova.cca(capscale.gen)
print(perm_anova.gen)

#repeat the same for Col-0 data

### calculate simpson diversity index

trans_otu<-t(otu_table)

#simpson diversity index
simpson<-diversity(trans_otu, "simpson")
simpson<-as.data.frame(simpson)
trans_otu_simpson <- cbind(design,simpson)

dfpwersimpson <- trans_otu_simpson %>% filter(compartment == "root" &
                                              genotype == "pWER::FLS2-GFP" &
                                              !biological.replicate == "c")

dfpwersimpson$cond<-as.factor(dfpwersimpson$cond)
levels(dfpwersimpson$cond)

ggplot(dfpwersimpson, aes(x=cond, y=simpson, fill=flg22)) +
  geom_boxplot(alpha=1, outlier.size=0, size=1.0, width=.8) +
  geom_jitter(position=position_jitterdodge(0.25), size=2, alpha=0.7) +
  scale_x_discrete(limits=c("NS FALSE pWER::FLS2-GFP", "NS TRUE pWER::FLS2-GFP", "S FALSE pWER::FLS2-GFP", "S TRUE pWER::FLS2-GFP","NS+S FALSE pWER::FLS2-GFP", "NS+S TRUE pWER::FLS2-GFP")) +
  scale_fill_manual(values=c("white","grey")) +
  theme_bw() +
  main_theme+
  labs(y="Simpson Index")+
  ggtitle("Simpson Index of pWER::FLS2-GFP plants")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggsave(paste(figures.dir, "Simpson_pwer.pdf", sep=""),  height=8, width=7)
#  ggsave(paste(figures.dir, "simpson_pwer.tiff", sep=""),  height=8, width=7)

# ANOVA
dfpwersimpson$SynCom_flg <- paste(dfpwersimpson$SynCom, dfpwersimpson$flg22)
dfpwersimpson$SynCom_flg<- as.factor(dfpwersimpson$SynCom_flg)
dfpwersimpson$biological.replicate <- design$biological.replicate[match(dfpwersimpson$SampleID, design$SampleID)]
dfpwersimpson$technical.replicate <- design$technical.replicate[match(dfpwersimpson$SampleID, design$SampleID)]

aov <-aov(simpson~ SynCom_flg+biological.replicate*technical.replicate, data=dfpwersimpson)
summary(aov)

tuk <- glht(aov, linfct = mcp(SynCom_flg = "Tukey"))
confint(tuk, level=0.95)
summary(tuk)
