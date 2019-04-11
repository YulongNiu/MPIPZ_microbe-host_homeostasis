# RNA-Seq data from Ka-Wai Ma #

<!-- content start -->

**Table of Contents**

- [1. Introduction](#1-introduction)
- [2. Tasks](#2-tasks)
    - [2.1 Alignment](#21-alignment)
    - [2.2 Quantification](#22-quantification)
- [3. Progress](#3-progress)
    - [3.1 pre-processing](#31-pre-processing)
    - [3.2 Alignment](#21-alignment)
- [4. Results](#4-results)
    - [4.1 Contamination](#41-contamination)
    - [4.2 Aligned reads](#42-aligned-reads)
- [5. Reference](#5-reference)
    
<!-- content end -->

## 1. Introduction ##

The *Arabidopsis thaliana* Col-0 treated with flg22, SynCom33, and SynCom35.

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

* Align reads to 

### 2.2 Quantification ###

* Quantification transcripts.

* Differentially expressed genes.

