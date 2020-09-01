# Guidance of RNA-Seq and amplicon sequencing analysis for the "microbe-host homeostasis" project

<!-- content start -->

- [1. RNA-Seq analysis](#1-rna-seq-analysis)
    - [1.1 Raw data](#11-raw-data)
    - [1.2 Pre-procession](#12-pre-procession)
- [2. Amplicon sequencing analysis](#2-amplicon-sequencing-analysis)
    - [2.1 Alignment](#21-alignment)
- [References](#references)

<!-- content end -->


## 1. RNA-Seq analysis

### 1.1 Raw data

Raw RNA-Seq data-sets can be retrieved from **GSE157128**. Arabidopsis thaliana* plants were germinated with the synthetic communities (SynComs) in the presence or absence of 1 µM flg22 and incubated for 14 days before harvesting. For *pWER::FLS2-GFP* plants, single-end RNA-Seq experiments were conducted for 4 conditions (3 biological replicates for each condition). For Col-0 plants, pair-end RNA-Seq experiments were conducted for 10 conditions (4 biological replicates for each condition).

### 1.2 Pre-procession

Raw Illumina RNA-Seq reads were pre-processed using fastp (v0.19.10) [1] with default settings for paired-end (Col-0 experiment) or single-end reads (*pWER::FLS2-GFP* experiment). For single-end reads, low quality sequences from the head (8 bases) and tail (2 bases) were trimmed. 

Scripts include [script_trim.sh](RNA-Seq_scripts/script_trim.sh).

### 1.3 


## 2. Amplicon sequencing analysis

## Reference

1. Chen, S., Zhou, Y., Chen, Y. & Gu, J. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics 34, i884-i890, doi:10.1093/bioinformatics/bty560 (2018).

2. 


