date

RAW_PATH=/biodata/dep_psl/grp_rgo/metatranscriptomics/data/flg22

links=(https://websafe.mpipz.mpg.de/d/TlcortL5BV/4096_A_run582_TCAGACGA_S1_L001_R1_001.fastq.gz
       https://websafe.mpipz.mpg.de/d/BtCt5Hukz5/4096_B_run582_GGCTTGTA_S2_L001_R1_001.fastq.gz
       https://websafe.mpipz.mpg.de/d/UfKDLFvCPO/4096_C_run582_GAGTCTCC_S3_L001_R1_001.fastq.gz
       https://websafe.mpipz.mpg.de/d/mhr10lZECs/4096_D_run582_CTAGACCT_S4_L001_R1_001.fastq.gz
       https://websafe.mpipz.mpg.de/d/bq7taBB0if/4096_E_run582_CGCAGTAC_S5_L001_R1_001.fastq.gz
       https://websafe.mpipz.mpg.de/d/pMtrZjcARt/4096_F_run582_CATGGACC_S6_L001_R1_001.fastq.gz
       https://websafe.mpipz.mpg.de/d/8Y21P2y7rr/4096_G_run582_ATTCGTCG_S7_L001_R1_001.fastq.gz
       https://websafe.mpipz.mpg.de/d/hL5SpVJmD1/4096_H_run582_AGGACTGT_S8_L001_R1_001.fastq.gz
       https://websafe.mpipz.mpg.de/d/t9CbsQIfI9/4096_I_run582_ACTCAAGT_S9_L001_R1_001.fastq.gz
       https://websafe.mpipz.mpg.de/d/2gpbn6EM4u/4096_J_run582_TTGCCATA_S10_L001_R1_001.fastq.gz
       https://websafe.mpipz.mpg.de/d/vwWbKQ8FMY/4096_K_run582_TATCCTCG_S11_L001_R1_001.fastq.gz
       https://websafe.mpipz.mpg.de/d/crLPuq4qJp/4096_L_run582_GGACTATT_S12_L001_R1_001.fastq.gz
       https://websafe.mpipz.mpg.de/d/KXfsO56Knk/4096_A_run582_TCAGACGA_S1_L002_R1_001.fastq.gz
       https://websafe.mpipz.mpg.de/d/leZP5K2L36/4096_B_run582_GGCTTGTA_S2_L002_R1_001.fastq.gz
       https://websafe.mpipz.mpg.de/d/XoltnAcfmz/4096_C_run582_GAGTCTCC_S3_L002_R1_001.fastq.gz
       https://websafe.mpipz.mpg.de/d/Ytnn6YDK26/4096_D_run582_CTAGACCT_S4_L002_R1_001.fastq.gz
       https://websafe.mpipz.mpg.de/d/Fxk2DhQf6j/4096_E_run582_CGCAGTAC_S5_L002_R1_001.fastq.gz
       https://websafe.mpipz.mpg.de/d/swApgB1BiM/4096_F_run582_CATGGACC_S6_L002_R1_001.fastq.gz
       https://websafe.mpipz.mpg.de/d/SWL6qnIV0Y/4096_G_run582_ATTCGTCG_S7_L002_R1_001.fastq.gz
       https://websafe.mpipz.mpg.de/d/S8iyOCGSGX/4096_H_run582_AGGACTGT_S8_L002_R1_001.fastq.gz
       https://websafe.mpipz.mpg.de/d/VTdhEFjLg6/4096_I_run582_ACTCAAGT_S9_L002_R1_001.fastq.gz
       https://websafe.mpipz.mpg.de/d/czYusmyYJe/4096_J_run582_TTGCCATA_S10_L002_R1_001.fastq.gz
       https://websafe.mpipz.mpg.de/d/oUg40ym2FA/4096_K_run582_TATCCTCG_S11_L002_R1_001.fastq.gz
       https://websafe.mpipz.mpg.de/d/AxxXWr2Bxf/4096_L_run582_GGACTATT_S12_L002_R1_001.fastq.gz)

cd ${RAW_PATH}

for i in ${links[@]}; do
    wget ${i}
done

date






date

