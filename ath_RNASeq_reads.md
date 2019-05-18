# How many reads are enough for Arabidopsis RNA-Seq

<!-- content start -->

**Table of Contents**

- [1. Estimation](#1-Estimation)
    - [1.1 Exons](#11-exons)
    - [1.2 Transcripts](#12-transcripts)
- [2. Expriments](#2-Estimation)

<!-- content end -->

## 1. Estimation

### 1.1 Exons

Estimate coverage from exons. 

1. Retrieve exons from [Ensembl database](ftp://ftp.ensemblgenomes.org/pub/plants/release-42/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz) (Araport11). A total of 193,130 non-redundant exons with average length of 336.4bp.

```
  Coverage `Reads (million)`
1 0.1X                0.0433
2 0.5X                0.217 
3 1X                  0.433 
4 2X                  0.866 
5 3X                  1.30  
6 5X                  2.17  
7 10X                 4.33  
8 15X                 6.50  
9 20X                 8.66  
```


