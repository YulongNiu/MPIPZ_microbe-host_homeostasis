### Scripts originally created by Ruben Garrido-Oter and modified by Ka-Wai Ma
# The scipt can be used to generate the box plots. Change the name of the files and samples e.g. strain ID accordingly for all 3 experiments

rm(list=ls())

if (!require(tidyverse)) install.packages("tidyverse")
if (!require(multcomp)) install.packages("multcomp")
library(tidyverse)
library(multcomp)

# set directories
results.dir <- "../results/"
data.dir <- "../data/"
figures.dir <- "../figures/"

source("plotting_functions.R")
source("cpcoa.func.R")

#change the name of the file for different experiments accordingly
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

# this is only necessary if we do spiking with E. coli sample
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


#generate a boxplot without subsetting first
df <- melt(otu_table_norm)
colnames(df) <- c("strain", "SampleID", "RA")

#match sampleID with the corresponding genotype, flg22, compartment, SynCom and cond
df$genotype <- design$genotype[match(df$SampleID, design$SampleID)]
df$flg22 <- design$flg22[match(df$SampleID, design$SampleID)]
df$compartment <- design$compartment[match(df$SampleID, design$SampleID)]
df$SynCom <- design$SynCom[match(df$SampleID, design$SampleID)]
df$cond <- design$cond[match(df$SampleID, design$SampleID)]
df$syncom <- taxonomy$syncom[match(df$strain, taxonomy$strain)]
df$biological.replicate <- design$biological.replicate[match(df$SampleID, design$SampleID)]
df$technical.replicate <- design$technical.replicate[match(df$SampleID, design$SampleID)]

#generate a table for colors and shapes
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
  c_very_light_red,
  c_very_light_red,
  c_very_light_red,
  c_very_light_red,
  c_cyan3,
  c_cyan3,
  c_cyan3,
  c_cyan3,
  c_green,
  c_green,
  c_green,
  c_green
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

colors <- colors[colors$group %in% df$cond, ]
df$cond <- factor(df$cond, levels=colors$group)

shapes <- shapes[shapes$group %in% df$cond, ]
df$cond <- factor(df$cond, levels=shapes$group)

df$family <- taxonomy$Family[match(df$strain, taxonomy$strain)]
df$syncomsource <- taxonomy$syncom[match(df$strain, taxonomy$strain)]

df$syncomsource_flg22 <- paste(df$syncomsource, df$flg22)
df$syncomsource_flg22 <- as.factor(df$syncomsource_flg22)

df$strain_syncom <- design$SynCom[match(df$SampleID, design$SampleID)]

df$strain_flg <- paste(df$strain, df$flg22)
df$strain_flg<-as.factor(df$strain_flg)

# we will use all data point except for replicate c in experiment 1 due to
# potential contamination or PCR error

NSpwer <- df %>% filter(compartment == "root" & genotype == "pWER::FLS2-GFP" &
                        strain_syncom =="NS")
NSpwer <- NSpwer %>% filter(!biological.replicate =="c")

str(NSpwer)
NSpwer$strain<-as.factor(NSpwer$strain)
NSpwer$syncomsource_flg22<-as.factor(NSpwer$syncomsource_flg22)
levels (NSpwer$strain)

Spwer <- df %>% filter(compartment == "root" & genotype == "pWER::FLS2-GFP" &
                       strain_syncom =="S")
Spwer <- Spwer %>% filter(!biological.replicate =="c")

str(Spwer)
Spwer$strain<-as.factor(Spwer$strain)
Spwer$syncomsource_flg22<-as.factor(Spwer$syncomsource_flg22)
levels (Spwer$strain)

NSSpwer <- df %>% filter(compartment == "root" & genotype == "pWER::FLS2-GFP" &
                         strain_syncom =="NS+S") 
NSSpwer <- NSSpwer %>% filter(!biological.replicate =="c")

str(NSSpwer)
NSSpwer$strain<-as.factor(NSSpwer$strain)
NSSpwer$syncomsource_flg22<-as.factor(NSSpwer$syncomsource_flg22)
levels (NSSpwer$strain)

NScol <- df %>% filter(compartment == "root" & genotype == "Col-0" &
                       strain_syncom =="NS")
NScol <- NScol %>% filter(!biological.replicate =="c")

str(NScol)
NScol$strain<-as.factor(NScol$strain)
NScol$syncomsource_flg22<-as.factor(NScol$syncomsource_flg22)

Scol <- df %>% filter(compartment == "root" & genotype == "Col-0" & strain_syncom =="S")
Scol <- Scol %>% filter(!biological.replicate =="c")

str(Scol)
Scol$strain<-as.factor(Scol$strain)
Scol$syncomsource_flg22<-as.factor(Scol$syncomsource_flg22)

NSScol <- df %>% filter(compartment == "root" & genotype == "Col-0" &
                        strain_syncom =="NS+S")
NSScol <- NSScol %>% filter(!biological.replicate =="c")

str(NSScol)
NSScol$strain<-as.factor(NSScol$strain)
NSScol$syncomsource_flg22<-as.factor(NSScol$syncomsource_flg22)


#Set main theme
main_theme <- theme(axis.line.x=element_line(color="black"),
                    axis.line.y=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(colour="black", size=10),
                    axis.text.x=element_text(colour="black", size=14,
                                             hjust = 1,vjust = 1,angle=45),
                    #use none to remove legend, can choose top, right, left
                    legend.position="top",
                    legend.title=element_blank(),
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    text=element_text(family="sans", size=15))

ggplot(NSpwer, aes(x=strain, y=RA, color=SynCom, shape=cond, fill=flg22)) +
  geom_boxplot(alpha=1, outlier.size=0, size=1.0, width=.8) +
  geom_jitter(position=position_jitterdodge(0.25), size=3, alpha=0.7) +
  #stat_summary(fun.y=mean, geom="point", shape=3, size=6.0) +
  scale_x_discrete(limits=c("Root60", "Root1290", "Root231","Root434","Root9")) +
  scale_colour_manual(values = c("#FF696CB3")) +
  scale_shape_manual(values=shapes$shape) +
  scale_fill_manual(values=c("white","#ffb3b4")) +
  ylim(0,0.85)+
  theme_bw() +
  main_theme +
  labs(y="Relative abundance (%)")+
  ggtitle("pWER::FLS2-GFP NS") +
  guides(color = FALSE, shape = FALSE) +
#  ggsave(paste(figures.dir, "box_NSpwerroot.tiff", sep=""), width=4, height=7) +
  ggsave(paste(figures.dir, "box_NSpwerroot.pdf", sep=""), width=4, height=7)


ggplot(Spwer, aes(x=strain, y=RA, color=SynCom, shape=cond,fill=flg22)) +
  geom_boxplot(alpha=1, outlier.size=0, size=1.0, width=.8) +
  geom_jitter(position=position_jitterdodge(0.25), size=3, alpha=0.7) +
  scale_x_discrete(limits=c("Root553", "Root342", "Root1212", "Root16D2", "Root179")) +
  scale_colour_manual(values = c("#00CDCDB3")) +
  scale_shape_manual(values=shapes$shape) +
  scale_fill_manual(values=c("white","#00e6e6")) +
  ylim(0,0.85)+
  theme_bw() +
  main_theme +
  labs(y="Relative abundance (%)")+
  ggtitle("pWER::FLS2-GFP S") + 
  guides(color = FALSE, shape = FALSE) +
#  ggsave(paste(figures.dir, "box_Spwerroot.tiff", sep=""), width=4, height=7) +
  ggsave(paste(figures.dir, "box_Spwerroot.pdf", sep=""), width=4, height=7)

ggplot(NSSpwer, aes(x=strain, y=RA, color=syncomsource, shape=cond, fill=syncomsource_flg22)) +
  geom_boxplot(alpha=1, outlier.size=0, size=1.0, width=.8 ) +
  geom_jitter(position=position_jitterdodge(0.25), size=3, alpha=0.7) +
  scale_x_discrete(limits=c("Root60","Root553", "Root1290","Root342", "Root231","Root1212","Root434","Root16D2","Root9", "Root179")) +
  scale_shape_manual(values=shapes$shape) +
  scale_fill_manual(values=c("white","#ffb3b4","white","#00e6e6")) +
  ylim(0,0.85)+
  theme_bw() +
  main_theme +
  labs(y="Relative abundance (%)")+
  ggtitle("pWER::FLS2-GFP NSS") +
  guides(color = FALSE, shape = FALSE) +
#  ggsave(paste(figures.dir, "box_NSSpwerroot.tiff", sep=""), width=9, height=7) +
  ggsave(paste(figures.dir, "box_NSSpwerroot.pdf", sep=""), width=9, height=7)
  

ggplot(NScol, aes(x=strain, y=RA, color=SynCom, shape=cond, fill=flg22)) +
  geom_boxplot(alpha=1, outlier.size=0, size=1.0, width=.8) +
  geom_jitter(position=position_jitterdodge(0.25), size=3, alpha=0.7) +
  scale_x_discrete(limits=c("Root60", "Root1290", "Root231","Root434","Root9")) +
  scale_colour_manual(values = c("#FF696CB3")) +
  scale_shape_manual(values=shapes$shape) +
  scale_fill_manual(values=c("white","#ffb3b4")) +
  ylim(0,0.85)+
  theme_bw() +
  main_theme +
  labs(y="Relative abundance (%)")+
  ggtitle("Col-0 NS") +
  guides(color = FALSE, shape = FALSE) +
#  ggsave(paste(figures.dir, "box_NScolroot.tiff", sep=""), width=4, height=7) +
  ggsave(paste(figures.dir, "box_NScolroot.pdf", sep=""), width=4, height=7)


ggplot(Scol, aes(x=strain, y=RA, color=SynCom, shape=cond,fill=flg22)) +
  geom_boxplot(alpha=1, outlier.size=0, size=1.0, width=.8) +
  geom_jitter(position=position_jitterdodge(0.25), size=3, alpha=0.7) +
  scale_x_discrete(limits=c("Root553", "Root342", "Root1212", "Root16D2", "Root179")) +
  scale_colour_manual(values = c("#00CDCDB3")) +
  scale_shape_manual(values=shapes$shape) +
  scale_fill_manual(values=c("white","#00e6e6")) +
  ylim(0,0.85)+
  theme_bw() +
  main_theme +
  labs(y="Relative abundance (%)")+
  ggtitle("Col-0 S") + 
  guides(color = FALSE, shape = FALSE) +
#  ggsave(paste(figures.dir, "box_Scolroot.tiff", sep=""), width=4, height=7) +
  ggsave(paste(figures.dir, "box_Scolroot.pdf", sep=""), width=4, height=7)

ggplot(NSScol, aes(x=strain, y=RA, color=syncomsource, shape=cond, fill=syncomsource_flg22)) +
  geom_boxplot(alpha=1, outlier.size=0, size=1.0, width=.8 ) +
  geom_jitter(position=position_jitterdodge(0.25), size=3, alpha=0.7) +
  scale_x_discrete(limits=c("Root60","Root553", "Root1290","Root342", "Root231","Root1212","Root434","Root16D2","Root9", "Root179")) +
  scale_shape_manual(values=shapes$shape) +
  scale_fill_manual(values=c("white","#ffb3b4","white","#00e6e6")) +
  ylim(0,0.85)+
  theme_bw() +
  main_theme +
  labs(y="Relative abundance (%)")+
  ggtitle("Col-0 NSS") +
  guides(color = FALSE, shape = FALSE) 
#  ggsave(paste(figures.dir, "box_NSScolroot.tiff", sep=""), width=9, height=7) +
  ggsave(paste(figures.dir, "box_NSScolroot.pdf", sep=""), width=9, height=7)

###compute statistics for different treatements by ANOVA

aovNSpwer <-aov(RA ~ strain_flg+biological.replicate*technical.replicate, data=NSpwer)
aovSpwer <-aov(RA ~ strain_flg+biological.replicate*technical.replicate, data=Spwer)
aovNSSpwer <-aov(RA ~ strain_flg+biological.replicate*technical.replicate, data=NSSpwer)
aovNScol <-aov(RA ~ strain_flg+biological.replicate*technical.replicate, data=NScol)
aovScol <-aov(RA ~ strain_flg+biological.replicate*technical.replicate, data=Scol)
aovNSScol <-aov(RA ~ strain_flg+biological.replicate*technical.replicate, data=NSScol)

#only one example was shown
tuk <- glht(aovNSpwer, linfct = mcp(strain_flg = "Tukey"))
confint(tuk, level=0.95)
summary(tuk)
