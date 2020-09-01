# Guidance of RNA-Seq and amplicon sequencing analysis for the "microbe-host homeostasis" project

<!-- content start -->

- [1. RNA-Seq analysis](#1-rna-seq-analysis)
    - [1.1 Raw data](#11-raw-data)
    - [1.2 Pre-procession](#12-pre-procession)
    - [1.3 Pseudo-alignment](#13-pseudo-alignment)
    - [1.4 DEGs analysis](#14-degs-analysis)
    - [1.5 Clustering](#15-clustering)
    - [1.6 GO enrichment](#16-go-enrichment)
- [2. Amplicon sequencing analysis](#2-amplicon-sequencing-analysis)
    - [2.1 Alignment](#21-alignment)
- [References](#references)

<!-- content end -->


## 1. RNA-Seq analysis

### 1.1 Raw data

Raw RNA-Seq data-sets can be retrieved from **GSE157128**. *Arabidopsis thaliana* plants were germinated with the synthetic communities (SynComs) in the presence or absence of 1 µM flg22 and incubated for 14 days before harvesting. For *pWER::FLS2-GFP* plants, single-end RNA-Seq experiments were conducted for 4 conditions (3 biological replicates for each condition). For Col-0 plants, pair-end RNA-Seq experiments were conducted for 10 conditions (4 biological replicates for each condition).

### 1.2 Pre-procession

Raw Illumina RNA-Seq reads were pre-processed using fastp (v0.19.10) [1] with default settings for paired-end (Col-0 experiment) or single-end reads (*pWER::FLS2-GFP* experiment). For single-end reads, low quality sequences from the head (8 bases) and tail (2 bases) were trimmed. 

Scripts:

[script_trim.sh](RNA-Seq_scripts/script_trim.sh).

### 1.3 Pseudo-alignment

High quality reads were pseudo-aligned to TAIR 10 *Arabidopsis thaliana* transcriptome reference (Ensembl) [2] using kallisto (v0.46.1) [3]. 

Scripts:

* [script_alignment.sh](RNA-Seq_scripts/script_alignment.sh).

Results: 

* pseudo-alignment count tables are [align_data](align_data) for *pWER::FLS2-GFP* and [align_data_1stadd](align_data_1stadd) for Col-0.

### 1.4 DEGs analysis

After removal of low abundant transcripts that were absent in at least two replicates under each condition, count data were imported using the tximport [4] package.

Differential expression analyses were performed using the DESeq2 [5] package. Firstly, raw counts were normalized with respect to the library size (*rlog* function) and log2 transformed. We tested for sample batch effects by surrogate variable (SV) analysis using the sva [6] package. Significant SVs were automatically detected and integrated into the model for differential analyses. Principal component analysis (*prcomp* function) based on whole transcripts were performed and plotted to visualize the cluster and variance of biological replicates under each condition. Abundance of Arabidopsis latent virus-1 reads did not correlate with sample variances and therefore removed from downstream analyses. Pair-wise comparisons were designed as: (1) flg22 treatment effect only, (2) non- suppressive and suppressive SynCom effect only, (3) flg22 treatment plus SynCom effects, (4) living vs. heat-killed bacteria. Transcripts with fold-changes > 1.5 and adjusted p-value for multiple comparisons (Benjamini–Hochberg method) equal to or below 0.05 were considered significant.

The log2 scaled counts were normalized by the identified SVs using the limma [7] package (*removeBatchEffect* function), and transformed as median-centered z-score by transcripts (scaled counts, *scale* function).

scripts: 

* [script_DEG.R](RNA-Seq_scripts/script_DEG.R) for *pWER::FLS2-GFP*.

* [script_DEG_1stadd.R](RNA-Seq_scripts/script_DEG_1stadd.R) for Col-0.

Results:

* DEGs tables are [eachGroup_vs_Mock_k.csv](results/removeZero/eachGroup_vs_Mock_k.csv) for *pWER::FLS2-GFP* [eachGroup_vs_Mock_k_1stadd.csv](results/removeZero/eachGroup_vs_Mock_k_1stadd.csv) for Col-0.

* Scales counts are [degres_condi_Mock.RData](results/removeZero/degres_condi_Mock.RData) for *pWER::FLS2-GFP* and [degres_condi_Mock_1stadd.RData](results/removeZero/degres_condi_Mock_1stadd.RData) for Col-0.

* PCA plots are [PCA_sva.pdf](results/removeZero/PCA_sva.pdf) for *pWER::FLS2-GFP* and [PCA_1stadd_sva.pdf](results/removeZero/PCA_1stadd_sva.pdf) for Col-0.

### 1.5 Clustering

The z-scores were used to conduct *k*-means clustering for all transcripts. The cluster number (*k* = 10) was determined by sum of squared error and Akaike information criterion. Next, confirmed transcripts with similar expression patterns were grouped in the same cluster. Differentially expressed transcripts (3,718 in *pWER::FLS2-GFP* and 4,450 in Col-0 experiments) and cluster results were visualized using heatmaps generated by ComplexHeatmap [8] package.

Scripts:

* [script_cluster.R](RNA-Seq_scripts/script_cluster.R) for *pWER::FLS2-GFP* *k*-means clustering.

* [script_cluster_1stadd.R](RNA-Seq_scripts/script_cluster_1stadd.R) for Col-0 *k*-means clustering.

* [script_replot_heatmap.R](RNA-Seq_scripts/script_replot_heatmap.R) for *pWER::FLS2-GFP* heatmap plot.

* [script_replot_heatmap_1stadd.R](RNA-Seq_scripts/script_replot_heatmap_1stadd.R) for Col-0 heatmap plot.

Results:

* Heatmaps are [kmeans10_heatmap_sig_DEG2.pdf](results/removeZero/kmeans10_heatmap_sig_DEG2.pdf) for *pWER::FLS2-GFP* and [kmeans10_heatmap_1stadd_sig_DEG2.pdf](results/removeZero/kmeans10_heatmap_1stadd_sig_DEG2.pdf) for Col-0.

### 1.6 GO enrichment

Gene ontology (GO) enrichment for each cluster using the whole Arabidopsis transcriptome as background were performed with the goseq [9] package with the consideration of transcripts length. GO annotations were retrieved from the Gene Ontology Consortium [10-11] (September 2019). Significantly changed biological process GO terms (adjusted p-value < 0.05) were visualized in dot plots using the clusterProfiler [12] package.

Scripts:

* [script_path.R](RNA-Seq_scripts/script_path.R) for *pWER::FLS2-GFP* GO enrichment.

* [script_path_1stadd.R](RNA-Seq_scripts/script_path_1stadd.R) for Col-0 GO enrichment.

Results:

* [kmeans10_cp_BP_dotplot_10.pdf](results/removeZero/geneset/clusterbc/kmeans10_cp_BP_dotplot_10.pdf) for *pWER::FLS2-GFP* and [kmeans10_1stadd_cp_BP_dotplot_10.pdf](results/removeZero/geneset_1stadd/clusterbc/kmeans10_1stadd_cp_BP_dotplot_10.pdf) GO dot plots.

## 2. Amplicon sequencing analysis

## References

1. Chen, S., Zhou, Y., Chen, Y. & Gu, J. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics 34, i884-i890, doi:10.1093/bioinformatics/bty560 (2018).

2. Yates, A. D. et al. Ensembl 2020. Nucleic Acids Res 48, D682-D688, doi:10.1093/nar/gkz966 (2020). 

3. Bray, N. L., Pimentel, H., Melsted, P. & Pachter, L. Near-optimal probabilistic RNA-seq quantification. Nat Biotechnol 34, 525-527, doi:10.1038/nbt.3519 (2016).

4. Soneson, C., Love, M. I. & Robinson, M. D. Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. F1000Res 4, 1521, doi:10.12688/f1000research.7563.2 (2015).

5. Love, M. I., Huber, W. & Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol 15, 550, doi:10.1186/s13059-014-0550-8 (2014). 

6. Leek, J. T., Johnson, W. E., Parker, H. S., Jaffe, A. E. & Storey, J. D. The sva package for removing batch effects and other unwanted variation in high-throughput experiments. Bioinformatics 28, 882-883, doi:10.1093/bioinformatics/bts034 (2012).

7. Ritchie, M. E. et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Res 43, e47, doi:10.1093/nar/gkv007 (2015).

8. Gu, Z., Eils, R. & Schlesner, M. Complex heatmaps reveal patterns and correlations in multidimensional genomic data. Bioinformatics 32, 2847-2849, doi:10.1093/bioinformatics/btw313 (2016).

9. Young, M. D., Wakefield, M. J., Smyth, G. K. & Oshlack, A. Gene ontology analysis for RNA- seq: accounting for selection bias. Genome Biol 11, R14, doi:10.1186/gb-2010-11-2-r14 (2010). 

10. Ashburner, M. et al. Gene ontology: tool for the unification of biology. The Gene Ontology Consortium. Nat Genet 25, 25-29, doi:10.1038/75556 (2000). 

11. The Gene Ontology, C. The Gene Ontology Resource: 20 years and still GOing strong. Nucleic Acids Res 47, D330-D338, doi:10.1093/nar/gky1055 (2019). 

12. Yu, G., Wang, L. G., Han, Y. & He, Q. Y. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS 16, 284-287, doi:10.1089/omi.2011.0118 (2012).



