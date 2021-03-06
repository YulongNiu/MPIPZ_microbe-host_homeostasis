# RNA-Seq data from Ka-Wai Ma

<!-- content start -->

- [1. pWER agar plate](#1-pwer-agar-plate)
    - [1.1 Alignment](#11-alignment)
    - [1.2 DEGs](#12-degs)
    - [1.3 Gene-set](#13-gene-set)
- [2. Col0 agar plate](#2-col0-agar-plate)
    - [2.1 Alignment](#21-alignment)
    - [2.2 DEGs](#22-degs)
    - [2.3 Gene-set](#23-gene-set)
- [3. Col0 flow pot](#3-col0-flow-pot)
    - [3.1 Alignment](#31-alignment)
    - [3.2 DEGs](#32-degs)
    - [3.3 Gene-set](#33-gene-set)
- [4. Merge](#4-merge)
    - [4.1 WER Col0 agar plate](#41-wer-col0-agar-plate)
    - [4.2 Col0 agar flowpot](#42-col0-agar-flowpot)
    - [4.3 Col0 agar flowpot iron](#43-col0-agar-flowpot-iron)
    - [4.4 Col0 agar flowpot flg22](#44-col0-agar-flowpot-flg22)
- [References](#references)

<!-- content end -->

## 1. pWER agar plate

### 1.1 Alignment

| sample           | rawfq    | trimfq   | H_ath | K_ath |
|------------------|----------|----------|-------|-------|
| Mock_1           | 16978228 | 16881903 | 0.977 | 0.962 |
| Mock_2           | 18851317 | 18723047 | 0.979 | 0.966 |
| Mock_3           | 18313224 | 18218407 | 0.98  | 0.966 |
| Flg22_1          | 17512255 | 17422675 | 0.971 | 0.954 |
| Flg22_2          | 18047704 | 17950283 | 0.975 | 0.959 |
| Flg22_3          | 17922052 | 17826125 | 0.974 | 0.957 |
| SynCom33_Flg22_1 | 17785611 | 17687693 | 0.974 | 0.955 |
| SynCom33_Flg22_2 | 16929833 | 16832003 | 0.975 | 0.96  |
| SynCom33_Flg22_3 | 20093432 | 19984770 | 0.974 | 0.956 |
| SynCom35_Flg22_1 | 17361283 | 17267286 | 0.965 | 0.952 |
| SynCom35_Flg22_2 | 16914038 | 16825207 | 0.979 | 0.966 |
| SynCom35_Flg22_3 | 21203518 | 21099083 | 0.972 | 0.957 |

### 1.2 DEGs

Remove `0|0|0` or `0|0|x`, then `29163` transcripts were kept.

* PCA plot

PCA plot of raw data

![PCA_raw](results/removeZero/PCA_raw.jpg)

PCA plot of SVA corrected data (3 hidden variables)

![PCA_sva](results/removeZero/PCA_sva.jpg)

* Cluster

![kmeans10](results/removeZero/kmeans10.jpg)

![kmeans10_trait](results/removeZero/kmeans10_trait.jpg)

* DEGs groups

```r
cond <- list(c('Mock_Flg22', 'Mock'),
             c('SynCom33_Flg22', 'Mock'),
             c('SynCom35_Flg22', 'Mock'),
             c('SynCom33_Flg22', 'Mock_Flg22'),
             c('SynCom35_Flg22', 'Mock_Flg22'),
             c('SynCom33_Flg22', 'SynCom35_Flg22'))
```

DEGs threshold: `|log2FC| > log2(1.5)` and `adj-pvalue < 0.05`

[DEGs table](results/removeZero/kmeans10_sig.csv)

* Heatmap

![kmeans10_heatmap_sig](results/removeZero/kmeans10_heatmap_sig.jpg)

![kmeans10_heatmap_sig2](results/removeZero/kmeans10_heatmap_sig2.jpg)

![kmeans10_heatmap_sig_DEG.jpg](results/removeZero/kmeans10_heatmap_sig_DEG.jpg)

![kmeans10_heatmap_sig_DEG.jpg2](results/removeZero/kmeans10_heatmap_sig_DEG2.jpg)

### 1.3 Gene-set

[GO analysis](results/removeZero/geneset)

* dotplot 

![kmeans10_cp_BP_dotplot](results/removeZero/geneset/clusterbc/kmeans10_cp_BP_dotplot.jpg)

![kmeans10_cp_BP_dotplot_10](results/removeZero/geneset/clusterbc/kmeans10_cp_BP_dotplot_10.jpg)

* Network (top 5 significant GO terms)

![kmeans10_cp_BP_network](results/removeZero/geneset/clusterbc/kmeans10_cp_BP_network.jpg)

* Optimized network (top 5 significant GO terms)

![kmeans10_cp_BP_cytoscape_5](results/removeZero/geneset/clusterbc/kmeans10_cp_BP_cytoscape_5.jpeg)

* Expand GO terms (top 10 significant GO terms)

![kmeans10_cp_BP_network_10](results/removeZero/geneset/clusterbc/kmeans10_cp_BP_network_10.jpg)

* Optimized network (top 10 significant GO terms)

![kmeans10_cp_BP_cytoscape_10](results/removeZero/geneset/clusterbc/kmeans10_cp_BP_cytoscape_10.jpeg)

## 2. Col0 agar plate

### 2.1 Alignment

| sample             | rawfq_R1 | rawfq_R2 | trimfq   | H_ath | K_ath | Hvirus_ath | Kvirus_ath | 
|--------------------|----------|----------|----------|-------|-------|------------|------------| 
| HKSynCom33_1       | 13172195 | 13172195 | 12783484 | 0.465 | 0.458 | 0.975      | 0.981      | 
| HKSynCom33_2       | 10767044 | 10767044 | 10287182 | 0.652 | 0.643 | 0.976      | 0.975      | 
| HKSynCom33_3       | 10878033 | 10878033 | 10560115 | 0.58  | 0.572 | 0.974      | 0.976      | 
| HKSynCom33_4       | 8651910  | 8651910  | 8454418  | 0.83  | 0.816 | 0.98       | 0.969      | 
| HKSynCom33_Flg22_1 | 6075012  | 6075012  | 5922387  | 0.887 | 0.871 | 0.978      | 0.964      | 
| HKSynCom33_Flg22_2 | 11499848 | 11499848 | 11142660 | 0.552 | 0.542 | 0.975      | 0.975      | 
| HKSynCom33_Flg22_3 | 11608298 | 11608298 | 11242512 | 0.718 | 0.71  | 0.976      | 0.975      | 
| HKSynCom33_Flg22_4 | 8763571  | 8763571  | 8519488  | 0.624 | 0.615 | 0.979      | 0.976      | 
| HKSynCom35_1       | 10084160 | 10084160 | 9826277  | 0.809 | 0.796 | 0.976      | 0.966      | 
| HKSynCom35_2       | 6297490  | 6297490  | 6170945  | 0.912 | 0.9   | 0.98       | 0.969      | 
| HKSynCom35_3       | 9915057  | 9915057  | 9603304  | 0.727 | 0.722 | 0.977      | 0.977      | 
| HKSynCom35_4       | 12587610 | 12587610 | 12091544 | 0.563 | 0.55  | 0.971      | 0.973      | 
| HKSynCom35_Flg22_1 | 6873637  | 6873637  | 6588694  | 0.977 | 0.963 | 0.977      | 0.963      | 
| HKSynCom35_Flg22_2 | 6461144  | 6461144  | 6307990  | 0.814 | 0.803 | 0.979      | 0.971      | 
| HKSynCom35_Flg22_3 | 8719346  | 8719346  | 8516934  | 0.6   | 0.595 | 0.977      | 0.98       | 
| HKSynCom35_Flg22_4 | 12973482 | 12973482 | 12618526 | 0.538 | 0.526 | 0.977      | 0.975      | 
| Mock_1             | 13036385 | 13036385 | 12732898 | 0.559 | 0.549 | 0.979      | 0.978      | 
| Mock_2             | 11495909 | 11495909 | 11144820 | 0.752 | 0.745 | 0.976      | 0.973      | 
| Mock_3             | 14230460 | 14230460 | 13715287 | 0.51  | 0.503 | 0.974      | 0.978      | 
| Mock_4             | 12354390 | 12354390 | 11513757 | 0.592 | 0.584 | 0.922      | 0.922      | 
| Mock_Flg22_1       | 11018794 | 11018794 | 10700223 | 0.703 | 0.693 | 0.978      | 0.973      | 
| Mock_Flg22_2       | 11452528 | 11452528 | 11123514 | 0.644 | 0.637 | 0.975      | 0.975      | 
| Mock_Flg22_3       | 13355598 | 13355598 | 12868862 | 0.579 | 0.571 | 0.974      | 0.977      | 
| Mock_Flg22_4       | 10720255 | 10720255 | 10352975 | 0.8   | 0.786 | 0.978      | 0.968      | 
| SynCom33_1         | 10665251 | 10665251 | 10238121 | 0.729 | 0.718 | 0.976      | 0.971      | 
| SynCom33_2         | 13777849 | 13777849 | 13220235 | 0.502 | 0.497 | 0.975      | 0.981      | 
| SynCom33_3         | 12414664 | 12414664 | 12026972 | 0.668 | 0.664 | 0.977      | 0.978      | 
| SynCom33_4         | 12574124 | 12574124 | 12030929 | 0.634 | 0.623 | 0.975      | 0.972      | 
| SynCom33_Flg22_1   | 11354523 | 11354523 | 11037970 | 0.67  | 0.661 | 0.976      | 0.974      | 
| SynCom33_Flg22_2   | 11098796 | 11098796 | 10723034 | 0.64  | 0.629 | 0.977      | 0.973      | 
| SynCom33_Flg22_3   | 12743113 | 12743113 | 12297414 | 0.536 | 0.531 | 0.878      | 0.881      | 
| SynCom33_Flg22_4   | 9369409  | 9369409  | 9138305  | 0.671 | 0.66  | 0.977      | 0.972      | 
| SynCom35_1         | 11683102 | 11683102 | 11384694 | 0.636 | 0.623 | 0.978      | 0.974      | 
| SynCom35_2         | 10254853 | 10254853 | 9911636  | 0.825 | 0.814 | 0.977      | 0.97       | 
| SynCom35_3         | 8767980  | 8767980  | 8025669  | 0.802 | 0.797 | 0.972      | 0.972      | 
| SynCom35_4         | 12844974 | 12844974 | 11942099 | 0.617 | 0.606 | 0.973      | 0.971      | 
| SynCom35_Flg22_1   | 11048570 | 11048570 | 10706105 | 0.614 | 0.603 | 0.977      | 0.974      | 
| SynCom35_Flg22_2   | 10988696 | 10988696 | 10515421 | 0.635 | 0.627 | 0.975      | 0.975      | 
| SynCom35_Flg22_3   | 8719287  | 8719287  | 8476158  | 0.707 | 0.698 | 0.977      | 0.974      | 
| SynCom35_Flg22_4   | 11436243 | 11436243 | 11120876 | 0.689 | 0.672 | 0.979      | 0.967      | 

### 2.2 DEGs

Remove `0|0|0|0` or `0|0|0|x`, then `27592` transcripts were kept.

[DEGs table](results/removeZero/kmeans10_1stadd_sig.csv)

* PCA plot

PCA plot of raw data

![PCA_1stadd_raw](results/removeZero/PCA_1stadd_raw.jpg)

PCA plot of SVA corrected data (3 hidden variables)

![PCA_1stadd_sva](results/removeZero/PCA_1stadd_sva.jpg)

* Cluster

![kmeans10_1stadd](results/removeZero/kmeans10_1stadd.jpg)

![kmeans10_trait_1stadd](results/removeZero/kmeans10_trait_1stadd.jpg)

* DEGs groups

```r
## flg22 treatment effect
cond1 <- list(c('Mock_Flg22', 'Mock'),
             c('HKSynCom33_Flg22', 'HKSynCom33'),
             c('HKSynCom35_Flg22', 'HKSynCom35'),
             c('SynCom33_Flg22', 'SynCom33'),
             c('SynCom35_Flg22', 'SynCom35'))

## bacteria effect
cond2 <- list(c('HKSynCom33', 'Mock'),
              c('HKSynCom35', 'Mock'),
              c('SynCom33', 'Mock'),
              c('SynCom35', 'Mock'),
              c('SynCom33', 'SynCom35'),
              c('HKSynCom33', 'HKSynCom35'))

## bacteria x flg22 effect
cond3 <- list(c('HKSynCom33_Flg22', 'Mock'),
              c('HKSynCom35_Flg22', 'Mock'),
              c('SynCom33_Flg22', 'Mock'),
              c('SynCom35_Flg22', 'Mock'),
              c('SynCom33_Flg22', 'SynCom35_Flg22'),
              c('HKSynCom33_Flg22', 'HKSynCom35_Flg22'),
              c('HKSynCom33_Flg22', 'Mock_Flg22'),
              c('HKSynCom35_Flg22', 'Mock_Flg22'),
              c('SynCom33_Flg22', 'Mock_Flg22'),
              c('SynCom35_Flg22', 'Mock_Flg22'))

## heat kill effect
cond4 <- list(c('SynCom33', 'HKSynCom33'),
              c('SynCom35', 'HKSynCom35'))
```

DEGs threshold: `|log2FC| > log2(1.5)` and `adj-pvalue < 0.05`

* Heatmap

![kmeans10_heatmap_1stadd_sig](results/removeZero/kmeans10_heatmap_1stadd_sig.jpg)

![kmeans10_heatmap_1stadd_sig2](results/removeZero/kmeans10_heatmap_1stadd_sig2.jpg)

![kmeans10_heatmap_1stadd_sig_DEG](results/removeZero/kmeans10_heatmap_1stadd_sig_DEG.jpg)

![kmeans10_heatmap_1stadd_sig_DEG](results/removeZero/kmeans10_heatmap_1stadd_sig_DEG2.jpg)

![kmeans10_heatmap_1stadd_sig_mergeMockHK_DEG](results/removeZero/kmeans10_heatmap_1stadd_sig_mergeMockHK_DEG.jpg)

![kmeans10_heatmap_1stadd_sig_mergeMockHK_DEG2](results/removeZero/kmeans10_heatmap_1stadd_sig_mergeMockHK_DEG2.jpg)

### 2.3 Gene-set

[GO analysis](results/removeZero/geneset_1stadd)

* dotplot 

![kmeans10_1stadd_cp_BP_dotplot](results/removeZero/geneset_1stadd/clusterbc/kmeans10_1stadd_cp_BP_dotplot.jpg)

* Network (top 5 significant GO terms)

![kmeans10_1stadd_cp_BP_network](results/removeZero/geneset_1stadd/clusterbc/kmeans10_1stadd_cp_BP_network.jpg)

* Optimized network (top 5 significant GO terms)

![kmeans10_1stadd_cp_BP_cytoscape_5](results/removeZero/geneset_1stadd/clusterbc/kmeans10_1stadd_cp_BP_cytoscape_5.jpeg)

* Expand GO terms (top 10 significant GO terms)

![kmeans10_1stadd_cp_BP_network_10](results/removeZero/geneset_1stadd/clusterbc/kmeans10_1stadd_cp_BP_network_10.jpg)

* Optimized network (top 10 significant GO terms)

![kmeans10_1stadd_cp_BP_cytoscape_10](results/removeZero/geneset_1stadd/clusterbc/kmeans10_1stadd_cp_BP_cytoscape_10.jpeg)

## 3. Col0 flow pot

### 3.1 Alignment

| sample     | rawfq_R1 | rawfq_R2 | trimfq   | H_ath | K_ath | Hvirus_ath | Kvirus_ath | 
|------------|----------|----------|----------|-------|-------|------------|------------| 
| Mock_1     | 8146907  | 8146907  | 7877992  | 0.743 | 0.733 | 0.967      | 0.963      | 
| Mock_2     | 18083448 | 18083448 | 17391986 | 0.245 | 0.241 | 0.96       | 0.98       | 
| Mock_3     | 12919098 | 12919098 | 12363313 | 0.468 | 0.465 | 0.962      | 0.974      | 
| SynCom33_1 | 13209898 | 13209898 | 12648097 | 0.406 | 0.401 | 0.959      | 0.97       | 
| SynCom33_2 | 15867598 | 15867598 | 15210307 | 0.277 | 0.274 | 0.968      | 0.986      | 
| SynCom33_3 | 12568370 | 12568370 | 12105344 | 0.437 | 0.435 | 0.968      | 0.98       | 
| SynCom35_1 | 10898141 | 10898141 | 10483412 | 0.652 | 0.649 | 0.969      | 0.975      | 
| SynCom35_2 | 16818587 | 16818587 | 16215582 | 0.341 | 0.338 | 0.965      | 0.98       | 
| SynCom35_3 | 11391744 | 11391744 | 10974092 | 0.57  | 0.566 | 0.967      | 0.973      | 

### 3.2 DEGs

Remove `0|0|0` or `0|0|x`, then `27105` transcripts were kept.

[DEGs table](results/removeZero/kmeans10_soil_sig.csv)

* PCA plot

PCA plot of raw data

![PCA_soil_raw](results/removeZero/PCA_soil_raw.jpg)

PCA plot of SVA corrected data (3 hidden variables)

![PCA_soil_sva](results/removeZero/PCA_soil_sva.jpg)

* Cluster

![kmeans10_soil](results/removeZero/kmeans10_soil.jpg)

![kmeans10_trait_soil](results/removeZero/kmeans10_trait_soil.jpg)

* DEGs groups

```r
cond <- list(c('SynCom33', 'Mock'),
             c('SynCom35', 'Mock'),
             c('SynCom33', 'SynCom35'))
```

* Heatmap

![kmeans10_heatmap_soil_sig](results/removeZero/kmeans10_heatmap_soil_sig.jpg)

![kmeans10_heatmap_soil_sig2](results/removeZero/kmeans10_heatmap_soil_sig2.jpg)

![kmeans10_heatmap_soil_sig_DEG](results/removeZero/kmeans10_heatmap_soil_sig_DEG.jpg)

![kmeans10_heatmap_soil_sig_DEG2](results/removeZero/kmeans10_heatmap_soil_sig_DEG2.jpg)

### 3.3 Gene-set

[GO analysis](results/removeZero/geneset_soil)

* dotplot 

![kmeans10_soil_cp_BP_dotplot_5](results/removeZero/geneset_soil/clusterbc/kmeans10_soil_cp_BP_dotplot_5.jpg)

![kmeans10_soil_cp_BP_dotplot_10](results/removeZero/geneset_soil/clusterbc/kmeans10_soil_cp_BP_dotplot_10.jpg)

* Network (top 5 significant GO terms)

![kmeans10_soil_cp_BP_network_5](results/removeZero/geneset_soil/clusterbc/kmeans10_soil_cp_BP_network_5.jpg)

* Optimized network (top 5 significant GO terms)

![kmeans10_soil_cp_BP_cytoscape_5](results/removeZero/geneset_soil/clusterbc/kmeans10_soil_cp_BP_cytoscape_5.jpeg)

* Expand GO terms (top 10 significant GO terms)

![kmeans10_soil_cp_BP_network_10](results/removeZero/geneset_soil/clusterbc/kmeans10_soil_cp_BP_network_10.jpg)

* Optimized network (top 10 significant GO terms)

![kmeans10_soil_cp_BP_cytoscape_10](results/removeZero/geneset_soil/clusterbc/kmeans10_soil_cp_BP_cytoscape_10.jpeg)

## 4. Merge

### 4.1 WER Col0 agar plate

![kmeans10_heatmap_WER_Col0](results/removeZero/kmeans10_heatmap_WER_Col0.jpg)

![kmeans10_heatmap_WER_Col02](results/removeZero/kmeans10_heatmap_WER_Col02.jpg)

### 4.2 Col0 agar flowpot

![kmeans10_heatmap_soil_agar](results/removeZero/kmeans10_heatmap_soil_agar.jpg)

![kmeans10_heatmap_soil_agar2](results/removeZero/kmeans10_heatmap_soil_agar2.jpg)

### 4.3 Col0 agar flowpot iron

![kmeans10_heatmap_WER_Col02_Iron](results/removeZero/kmeans10_heatmap_WER_Col02_Iron.jpg)

![kmeans10_heatmap_WER_Col02_Iron2.jpg](results/removeZero/kmeans10_heatmap_WER_Col02_Iron2.jpg)

### 4.4 Col0 agar flowpot flg22

![kmeans10_heatmap_WER_Col02_flg22](results/removeZero/kmeans10_heatmap_WER_Col02_flg22.jpg)

![kmeans10_heatmap_WER_Col02_flg222.jpg](results/removeZero/kmeans10_heatmap_WER_Col02_flg222.jpg)

## References

* Garrido-Oter R, Nakano RT, Dombrowski N, Ma KW; AgBiome Team, McHardy AC, Schulze-Lefert P. **Modular Traits of the Rhizobiales Root Microbiota and Their Evolutionary Relationship with Symbiotic Rhizobia** *Cell Host Microbe.* 2018;24(1):155-167.e5.
