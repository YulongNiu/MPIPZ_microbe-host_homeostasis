# RNA-Seq data from Ka-Wai Ma #

<!-- content start -->

**Table of Contents**

- [1. Introduction](#1-introduction)
- [2. Tasks](#2-tasks)
    - [2.1 Alignment](#21-alignment)
    - [2.2 Quantification](#22-quantification)
- [3. Progress](#3-progress)
    - [3.1 Pre-processing](#31-pre-processing)
    - [3.2 Alignment](#21-alignment)
- [4. Results](#4-results)
    - [4.1 Aligned reads](#41-aligned-reads)
    - [4.2 DEGs](#42-DEGs)
    - [4.3 Cluster](#43-Cluster)
        - [4.3.1 Hierarchical clustering](#431-hierarchical-clustering)
        - [4.3.2 k-means clustering](#432-k-means-clustering)
    - [4.4 Gene-set analysis](#44-gene-set-analysis)
- [5. 1stadd](#5-1stadd)
    - [5.1 Aligned reads](#51-aligned-reads)
    - [5.2 DEGs](#52-DEGs)
    - [5.3 Cluster](#53-cluster)
- [6. Comparison](#6-comparison)
- [References](#references)
    
<!-- content end -->

## 1. Introduction ##

The *Arabidopsis thaliana* Col-0 treated with flg22, SynCom33+flg22, and SynCom35+flg22.

> reference genomes path `/netscratch/dep_psl/grp_rgo/yniu/ref`

> RNA-seq data path `/biodata/dep_psl/grp_rgo/metatranscriptomics/data/flg22`

| sample number | genotype At    | organ | flg22 | bacteria  | biological replicate | Library title | Library number | Library state | Reads requested | Reads sequenced | 
|---------------|----------------|-------|-------|-----------|----------------------|---------------|----------------|---------------|-----------------|-----------------| 
| a             | pWER::FLS2-GFP | roots | 0.00  | no        | 1                    | a.1           | 4016.A.2       | Draft         | 20,000,000      | 0.0             | 
| b             | pWER::FLS2-GFP | roots | 1.00  | no        | 1                    | b.1           | 4016.B.2       | Draft         | 20,000,000      | 0.0             | 
| c             | pWER::FLS2-GFP | roots | 1.00  | SynCom33  | 1                    | c.1           | 4016.C.2       | Draft         | 20,000,000      | 0.0             | 
| d             | pWER::FLS2-GFP | roots | 1.00  | SynCom35  | 1                    | d.1           | 4016.D.2       | Draft         | 20,000,000      | 0.0             | 
| e             | pWER::FLS2-GFP | roots | 0.00  | no        | 2                    | e.1           | 4016.E.2       | Draft         | 20,000,000      | 0.0             | 
| f             | pWER::FLS2-GFP | roots | 1.00  | no        | 2                    | f.1           | 4016.F.2       | Draft         | 20,000,000      | 0.0             | 
| g             | pWER::FLS2-GFP | roots | 1.00  | SynCom33  | 2                    | g.1           | 4016.G.2       | Draft         | 20,000,000      | 0.0             | 
| h             | pWER::FLS2-GFP | roots | 1.00  | SynCom35  | 2                    | h.1           | 4016.H.2       | Draft         | 20,000,000      | 0.0             | 
| i             | pWER::FLS2-GFP | roots | 0.00  | no        | 3                    | i.1           | 4016.I.2       | Draft         | 20,000,000      | 0.0             | 
| j             | pWER::FLS2-GFP | roots | 1.00  | no        | 3                    | j.1           | 4016.J.2       | Draft         | 20,000,000      | 0.0             | 
| k             | pWER::FLS2-GFP | roots | 1.00  | SynCom33  | 3                    | k.1           | 4016.K.2       | Draft         | 20,000,000      | 0.0             | 
| l             | pWER::FLS2-GFP | roots | 1.00  | SynCom35  | 3                    | l.1           | 4016.L.2       | Draft         | 20,000,000      | 0.0             | 

> analysis path `/netscratch/dep_psl/grp_rgo/yniu/KaWaiFlg22`

## 2. Tasks ##

### 2.1 Alignment ###

* Align reads to *Arabidopsis thaliana* Col-0 genomes (HISAT2) and cDNAs (Kallisto)

### 2.2 Quantification ###

* Quantification transcripts.

* Differentially expressed genes.

## 3. Progress ##

### 3.1 Pre-processing ###

* Trim by [fastp](https://github.com/OpenGene/fastp) with stringent parameters that remove the 8-mer in the head and 2-mer in the tail of each read.

### 3.2 Alignment ###

* Use Kallisto and HISAT2 to align reads.

## 4. Results ##

### 4.1 Aligned reads ###

* Alignment rates

| sample           | rawfq    | trimfq   | H_ath | K_ath | 
|------------------|----------|----------|-------|-------| 
| Mock_1           | 16978228 | 16881903 | 0.977 | 0.962 | 
| Mock_2           | 18851317 | 18723047 | 0.979 | 0.966 | 
| Mock_3           | 18313224 | 18218407 | 0.98  | 0.966 | 
| Flg22_1          | 17512255 | 17422675 | 0.971 | 0.954 | 
| Flg22_2          | 18047704 | 17950283 | 0.975 | 0.959 | 
| Flg22_3          | 17922052 | 17826125 | 0.974 | 0.957 | 
| Flg22_SynCom33_1 | 17785611 | 17687693 | 0.974 | 0.955 | 
| Flg22_SynCom33_2 | 16929833 | 16832003 | 0.975 | 0.96  | 
| Flg22_SynCom33_3 | 20093432 | 19984770 | 0.974 | 0.956 | 
| Flg22_SynCom35_1 | 17361283 | 17267286 | 0.965 | 0.952 | 
| Flg22_SynCom35_2 | 16914038 | 16825207 | 0.979 | 0.966 | 
| Flg22_SynCom35_3 | 21203518 | 21099083 | 0.972 | 0.957 | 

> "K" and "H" is for Kallisto and HISAT2 alignment.

* Kallisto alignment example

![K_alignment](results/Kallisto_alignment.png)

* HISAT2 alignment example

![H_alignment](results/HISAT2_alignment.png)

* PCA plot

![PCA](results/PCA.jpg)

> rlog (regularized log transformation) in PCA.

### 4.2 DEGs ###

```
                      condition
Mock_1                     Mock
Mock_2                     Mock
Mock_3                     Mock
Flg22_1                   Flg22
Flg22_2                   Flg22
Flg22_3                   Flg22
Flg22_SynCom33_1 Flg22_SynCom33
Flg22_SynCom33_2 Flg22_SynCom33
Flg22_SynCom33_3 Flg22_SynCom33
Flg22_SynCom35_1 Flg22_SynCom35
Flg22_SynCom35_2 Flg22_SynCom35
Flg22_SynCom35_3 Flg22_SynCom35
```

* 3 groups (flg22, flg22_SynCom33, and flg22_SynCom35) *vs.* Mock

[results/eachGroup_vs_Mock_k_full.csv](results/eachGroup_vs_Mock_k_full.csv) (remove genes with zero count) and [results/eachGroup_vs_Mock_k.csv](results/eachGroup_vs_Mock_k.csv) (remove gene with at least two zero count in one condition).

* 2 groups (flg22_SynCom33 and flg22_SynCom35) *vs.* flg22

[results/SynCom_vs_flg22_k_full.csv](results/SynCom_vs_flg22_k_full.csv) and [results/SynCom_vs_flg22_k.csv](results/SynCom_vs_flg22_k.csv).

* 1 groups flg22_SynCom35 *vs.* flg22_SynCom33

[results/SynCom35_vs_SynCom33_k_full.csv](SynCom35_vs_SynCom33_k_full.csv) and [results/SynCom35_vs_SynCom33_k.csv](results/SynCom35_vs_SynCom33_k.csv).

> logarithm transformation in columns `Mock_1` to `Flg22_SynCom35_3`.

### 4.3 Cluster ###

#### 4.3.1 Hierarchical clustering ####

Hierarchical clustering to find potential gene expression patterns in four conditions (`Mock`, `flg22`, `SynCom33`, and `SynCom35`). The read count of each transcript was scaled, then average value was used to represent each condition. Transcript with more than one zero counts in any of the four conditions were excluded (29161 transcripts left).

Distance was calculated as `1 - Pearson's Correlation Coefficient`

![genetree](results/genetree.jpg)

* Cut trees to get clusters

  height threshold `0.5`
  
  ![hieracluster_0d5](results/hieracluster_0d5.jpg)
  
  ```
  AT1G14550.1 AT2G30750.1 AT2G19190.1 
         21          26          20 
  ```

  height threshold `1.0`
  
  ![hieracluster_1d0](results/hieracluster_1d0.jpg)
  
  ```
  AT1G14550.1 AT2G30750.1 AT2G19190.1 
          9          14          14 
  ```
  
  height threshold `1.5`
  
  ```
  AT1G14550.1 AT2G30750.1 AT2G19190.1 
          7           7           7 
  ```
  
  ![hieracluster_1d5](results/hclust_1d5.jpg)
  
  ![hieracluster_gene_1d5](results/hclust_1d5_genes.jpg)
  
* correlation of traits (phenotype) and clusters
  
  ```
  flg22 SynCom33 SynCom35 bacteria rootlen
      0        0        0        0     5.6
      1        0        0        0     1.1
      1        1        0        1     1.3
      1        0        1        1     4.9
  ```
  
  ![hieracluster_1d5_trait](results/hclust_1d5_trait.jpg)
  
* heatmaps
  
  ![heatmap_merge](results/hclust_1d5_heatmap_merge.jpg)
  
  ![heatmap_group](results/hclust_1d5_heatmap_group.jpg)
  
  ![heatmap_raw](results/hclust_1d5_heatmap_raw.jpg)
  
  ![heatmap_scale](results/hclust_1d5_heatmap_scale.jpg)
  
  ![heatmap_logFC](results/hclust_1d5_heatmap_logFC.jpg)
  
  ![heatmap_sig](results/hclust_1d5_heatmap_sig.jpg)
  
#### 4.3.2 k-means clustering ####

Choose k (group number) by the sum of squared error (SSE) and Akaike information criterion (AIC). k number **10** was finally chose.

![kmeans_sse](results/kmeans_sse.jpg)

![kmeans_AIC](results/kmeans_AIC.jpg)

![kmeans_1d5](results/kmeans_10.jpg)
  
![kmeans_gene_1d5](results/kmeans_10_genes.jpg)

  ```
  AT1G14550.1 AT2G30750.1 AT2G19190.1 
          9           9           9 
  ```
* correlation of traits (phenotype) and clusters
  
  ```
  flg22 SynCom33 SynCom35 bacteria rootlen
      0        0        0        0     5.6
      1        0        0        0     1.1
      1        1        0        1     1.3
      1        0        1        1     4.9
  ```
  
  ![kmeans_1d5_trait](results/kmeans_10_trait.jpg)
  
* heatmaps
  
  ![heatmap_merge](results/kmeans_10_heatmap_merge.jpg)
  
  ![heatmap_group](results/kmeans_10_heatmap_group.jpg)
  
  ![heatmap_raw](results/kmeans_10_heatmap_raw.jpg)
  
  ![heatmap_scale](results/kmeans_10_heatmap_scale.jpg)
  
  ![heatmap_logFC](results/kmeans_10_heatmap_logFC.jpg)
  
  ![heatmap_sig](results/kmeans_10_heatmap_sig.jpg)
  
### 4.4 Gene-set analysis ###

Find significantly changed gene ontology, KEGG, and BioCyc gene-sets in each k-means cluster.

![means10_cluster9_BP](results/pathway/kmeans10_cluster9_BP.jpg)

![means10_cluster9_MF](results/pathway/kmeans10_cluster9_MF.jpg)

![means10_cluster9_CC](results/pathway/kmeans10_cluster9_CC.jpg)

![means10_cluster9_KEGG](results/pathway/kmeans10_cluster9_KEGG.jpg)

![means10_cluster9_BioCyc](results/pathway/kmeans10_cluster9_BioCyc.jpg)

## 5. 1stadd

### 5.1 aligned reads

| sample             | rawfq   | trimfq  | H_ath | K_ath | Hvirus_ath | Kvirus_ath | 
|--------------------|---------|---------|-------|-------|------------|------------| 
| HKSynCom33_1       | 6133228 | 6009879 | 0.469 | 0.461 | 0.977      | 0.982      | 
| HKSynCom33_2       | 6509281 | 6214550 | 0.654 | 0.645 | 0.977      | 0.975      | 
| HKSynCom33_3       | 6556180 | 6416355 | 0.583 | 0.574 | 0.976      | 0.976      | 
| HKSynCom33_4       | 8651910 | 8454418 | 0.831 | 0.816 | 0.98       | 0.969      | 
| HKSynCom33_Flg22_1 | 6075012 | 5922387 | 0.887 | 0.871 | 0.978      | 0.964      | 
| HKSynCom33_Flg22_2 | 6675351 | 6518231 | 0.555 | 0.545 | 0.976      | 0.975      | 
| HKSynCom33_Flg22_3 | 6949688 | 6786291 | 0.72  | 0.711 | 0.978      | 0.975      | 
| HKSynCom33_Flg22_4 | 8763571 | 8519488 | 0.624 | 0.615 | 0.979      | 0.976      | 
| HKSynCom35_1       | 5908989 | 5794978 | 0.812 | 0.798 | 0.978      | 0.967      | 
| HKSynCom35_2       | 6297490 | 6170945 | 0.912 | 0.9   | 0.98       | 0.969      | 
| HKSynCom35_3       | 5647852 | 5514059 | 0.73  | 0.724 | 0.978      | 0.977      | 
| HKSynCom35_4       | 5942014 | 5719804 | 0.566 | 0.553 | 0.972      | 0.974      | 
| HKSynCom35_Flg22_1 | 6873637 | 6588694 | 0.977 | 0.963 | 0.977      | 0.963      | 
| HKSynCom35_Flg22_2 | 6461144 | 6307990 | 0.814 | 0.803 | 0.979      | 0.971      | 
| HKSynCom35_Flg22_3 | 8719346 | 8516934 | 0.6   | 0.595 | 0.978      | 0.98       | 
| HKSynCom35_Flg22_4 | 8498860 | 8320269 | 0.54  | 0.528 | 0.979      | 0.975      | 
| Mock_1             | 8150057 | 8004876 | 0.561 | 0.551 | 0.981      | 0.979      | 
| Mock_2             | 7319673 | 7147533 | 0.754 | 0.746 | 0.977      | 0.974      | 
| Mock_3             | 5657217 | 5520195 | 0.515 | 0.506 | 0.976      | 0.978      | 
| Mock_4             | 8016274 | 7331269 | 0.626 | 0.618 | 0.973      | 0.975      | 
| Mock_Flg22_1       | 5753204 | 5629060 | 0.706 | 0.695 | 0.98       | 0.974      | 
| Mock_Flg22_2       | 6756972 | 6609999 | 0.647 | 0.639 | 0.977      | 0.975      | 
| Mock_Flg22_3       | 6095678 | 5941698 | 0.582 | 0.574 | 0.976      | 0.977      | 
| Mock_Flg22_4       | 6114701 | 5915059 | 0.801 | 0.788 | 0.978      | 0.968      | 
| SynCom33_1         | 6496075 | 6185955 | 0.73  | 0.72  | 0.975      | 0.971      | 
| SynCom33_2         | 6573180 | 6316721 | 0.505 | 0.5   | 0.975      | 0.981      | 
| SynCom33_3         | 8633739 | 8410426 | 0.671 | 0.665 | 0.979      | 0.978      | 
| SynCom33_4         | 5998570 | 5741137 | 0.636 | 0.626 | 0.975      | 0.973      | 
| SynCom33_Flg22_1   | 6049394 | 5933045 | 0.673 | 0.664 | 0.978      | 0.974      | 
| SynCom33_Flg22_2   | 6433414 | 6228966 | 0.643 | 0.632 | 0.977      | 0.974      | 
| SynCom33_Flg22_3   | 6083498 | 5921392 | 0.598 | 0.591 | 0.977      | 0.978      | 
| SynCom33_Flg22_4   | 9369409 | 9138305 | 0.671 | 0.66  | 0.977      | 0.972      | 
| SynCom35_1         | 7208495 | 7066347 | 0.638 | 0.625 | 0.98       | 0.974      | 
| SynCom35_2         | 6036827 | 5884003 | 0.828 | 0.816 | 0.978      | 0.97       | 
| SynCom35_3         | 8767980 | 8025669 | 0.802 | 0.797 | 0.972      | 0.972      | 
| SynCom35_4         | 7560214 | 6873333 | 0.618 | 0.608 | 0.971      | 0.972      | 
| SynCom35_Flg22_1   | 6224996 | 6077825 | 0.617 | 0.605 | 0.979      | 0.975      | 
| SynCom35_Flg22_2   | 6208664 | 6005766 | 0.639 | 0.629 | 0.976      | 0.976      | 
| SynCom35_Flg22_3   | 8719287 | 8476158 | 0.707 | 0.698 | 0.977      | 0.974      | 
| SynCom35_Flg22_4   | 6892972 | 6745962 | 0.691 | 0.673 | 0.98       | 0.968      | 

### 5.2 DEGs

![PCA_1stadd](results/PCA_1stadd.jpg)

![PCA_1stadd_3d](results/PCA_1stadd_3d.gif)

threshold 2

![PCA_1stadd](results/PCA_1stadd_thres2.jpg)

threshold 5

![PCA_1stadd](results/PCA_1stadd_thres5.jpg)

* 1st RNA-Seq

![kmeans_10_selectgene](results/kmeans_10_selectgene.jpg)

* 2nd RNA-Seq

![cluster10_selectgene_1stadd](results/cluster10_selectgene_1stadd.jpg)

### 5.3 Cluster

![kmeans_sse_1stadd](results/kmeans_sse_1stadd.jpg)

![kmeans_AIC_1stadd](results/kmeans_AIC_1stadd.jpg)

![cluster10_1stadd](results/cluster10_1stadd.jpg)

![cluster10_genes_1stadd](results/cluster10_genes_1stadd.jpg)

```
   Flg22 SynCom33 SynCom35 bacteria
1      0        0        0        0
2      1        0        0        0
3      0        1        0        0
4      1        1        0        0
5      0        1        0        1
6      1        1        0        1
7      0        0        1        0
8      1        0        1        0
9      0        0        1        1
10     1        0        1        1
```

![cluster10_trait_1stadd](results/cluster10_trait_1stadd.jpg)

![1stadd.jpg](results/1stadd.jpg)

* Mock HKSynCom33 HKSynCom35

![cluster10_g1_1stadd](results/cluster10_g1_1stadd.jpg)

* Mock Mock+flg22 HKSynCom33+flg22 HKSynCom35+flg22

![cluster10_g2_1stadd](results/cluster10_g2_1stadd.jpg)

* Mock SynCom33 SynCom35

![cluster10_g3_1stadd](results/cluster10_g3_1stadd.jpg)

* Mock Mock+flg22 SynCom33+flg22 SynCom35+flg22

![cluster10_g4_1stadd](results/cluster10_g4_1stadd.jpg)

## 6. Comparison 

* SynCom35 down DEGs

![veen_plot_35down](results/veen_plot_35down.jpg)

![]()

* SynCom35 up DEGs

![veen_plot_35up](results/veen_plot_35up.jpg)



## References

* Garrido-Oter R, Nakano RT, Dombrowski N, Ma KW; AgBiome Team, McHardy AC, Schulze-Lefert P. **Modular Traits of the Rhizobiales Root Microbiota and Their Evolutionary Relationship with Symbiotic Rhizobia** *Cell Host Microbe.* 2018;24(1):155-167.e5.

<!-- * [Clustering RNAseq data, making heatmaps, and tree cutting to identify gene modules.](https://2-bitbio.com/2017/04/clustering-rnaseq-data-making-heatmaps.html) -->






