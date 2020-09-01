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

### 1.3 Pseudo-alignment

High quality reads were pseudo-aligned to TAIR 10 *Arabidopsis thaliana* transcriptome reference (Ensembl) [2] using kallisto (v0.46.1) [3]. 

Scripts:

* [script_alignment.sh](RNA-Seq_scripts/script_alignment.sh).

Results: 

* pseudo-alignment count tables are [align_data](align_data) for *pWER::FLS2-GFP* and [align_data_1stadd](align_data_1stadd) for Col-0.

### 1.4 DEGs analysis

After removal of low abundant transcripts that were absent in at least two replicates under each condition, count data were imported using the tximport [4] package.

Differential expression analyses were performed using the DESeq2 [5] package. Firstly, raw counts were normalized with respect to the library size (rlog function) and log2 transformed. We tested for sample batch effects by surrogate variable (SV) analysis using the sva [6] package. Significant SVs were automatically detected and integrated into the model for differential analyses. Principal component analysis (prcomp function) based on whole transcripts were performed and plotted to visualize the cluster and variance of biological replicates under each condition. Abundance of Arabidopsis latent virus-1 reads did not correlate with sample variances and therefore removed from downstream analyses. Pair-wise comparisons were designed as: (1) flg22 treatment effect only, (2) non- suppressive and suppressive SynCom effect only, (3) flg22 treatment plus SynCom effects, (4) living vs. heat-killed bacteria. Transcripts with fold-changes > 1.5 and adjusted p-value for multiple comparisons (Benjamini–Hochberg method) equal to or below 0.05 were considered significant.

scripts: 

* [script_DEG.R](RNA-Seq_scripts/script_DEG.R) for *pWER::FLS2-GFP*.

* [script_DEG_1stadd.R](RNA-Seq_scripts/script_DEG_1stadd.R) for Col-0.

Results:

*  *pWER::FLS2-GFP* 

* for Col-0.

## 2. Amplicon sequencing analysis

## Reference

1. Chen, S., Zhou, Y., Chen, Y. & Gu, J. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics 34, i884-i890, doi:10.1093/bioinformatics/bty560 (2018).

2. Yates, A. D. et al. Ensembl 2020. Nucleic Acids Res 48, D682-D688, doi:10.1093/nar/gkz966 (2020). 

3. Bray, N. L., Pimentel, H., Melsted, P. & Pachter, L. Near-optimal probabilistic RNA-seq quantification. Nat Biotechnol 34, 525-527, doi:10.1038/nbt.3519 (2016).

4. Soneson, C., Love, M. I. & Robinson, M. D. Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. F1000Res 4, 1521, doi:10.12688/f1000research.7563.2 (2015).

5. Love, M. I., Huber, W. & Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol 15, 550, doi:10.1186/s13059-014-0550-8 (2014). 

6. Leek, J. T., Johnson, W. E., Parker, H. S., Jaffe, A. E. & Storey, J. D. The sva package for removing batch effects and other unwanted variation in high-throughput experiments. Bioinformatics 28, 882-883, doi:10.1093/bioinformatics/bts034 (2012).


