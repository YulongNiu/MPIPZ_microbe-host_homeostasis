
## originally by Yulong Niu
## yulong.niu@hotmail.com

date

REF_PATH=/netscratch/dep_psl/grp_rgo/yniu/ref/ath
CLEAN_PATH=/netscratch/dep_psl/grp_rgo/yniu/KaWaiFlg22/clean_data_1stadd
ALIGN_PATH=/netscratch/dep_psl/grp_rgo/yniu/KaWaiFlg22/align_data_1stadd_virus1
KALLISTO_PATH=/home/yniu/Biotools/kallisto_v0.46.1
HISAT2_PATH=/home/yniu/Biotools/hisat2-2.1.0
SAMTOOL_PATH=/opt/share/software/bin
RM_PATH=/bin/rm
MOVE_PATH=/bin/mv

CORENUM=45
SPECIES='ath'

cd ${REF_PATH}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~raw~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fq=($(ls ${CLEAN_PATH} | grep .fq.gz))
fqnames=($(echo "${fq[@]%%.*}" | tr ' ' '\n' | sort -u | tr '\n' ' '))

for i in ${fqnames[@]}; do

    echo "===================================="
    echo "Kallisto using ${SPECIES} cDNA for ${i}."
    ${KALLISTO_PATH}/kallisto quant \
                    --genomebam \
                    --gtf ${REF_PATH}/Arabidopsis_thaliana.TAIR10.42.gtf.gz \
                    --chromosomes ${REF_PATH}/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.fai \
                    -t ${CORENUM} \
                    -i ${REF_PATH}/ath.kindex \
                    -o ${ALIGN_PATH}/${i}_${SPECIES}_kallisto \
                    --single -l 200 -s 20\
                    ${CLEAN_PATH}/${i}.fq.gz

    ${MOVE_PATH} ${ALIGN_PATH}/${i}_${SPECIES}_kallisto/pseudoalignments.bam \
                 ${ALIGN_PATH}/${i}_${SPECIES}_kallisto/${i}_${SPECIES}_pseudoalignments.bam

    ${MOVE_PATH} ${ALIGN_PATH}/${i}_${SPECIES}_kallisto/pseudoalignments.bam.bai \
                 ${ALIGN_PATH}/${i}_${SPECIES}_kallisto/${i}_${SPECIES}_pseudoalignments.bam.bai

    echo "HISAT2 using ${SPECIES} genome for ${i}."
    ${HISAT2_PATH}/hisat2 -p ${CORENUM} \
                  -x ${REF_PATH}/athht2index/genome \
                  -U ${CLEAN_PATH}/${i}.fq.gz \
                  -S ${ALIGN_PATH}/${i}_${SPECIES}_hisat2.sam

    ${SAMTOOL_PATH}/samtools sort \
                   -@ ${CORENUM} \
                   -o ${ALIGN_PATH}/${i}_${SPECIES}_hisat2.bam \
                   ${ALIGN_PATH}/${i}_${SPECIES}_hisat2.sam

    ${SAMTOOL_PATH}/samtools index ${ALIGN_PATH}/${i}_${SPECIES}_hisat2.bam

    ${RM_PATH} ${ALIGN_PATH}/${i}_${SPECIES}_hisat2.sam
    echo "====================================="

done
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~1stadd~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fq=($(ls ${CLEAN_PATH} | grep .fq.gz))
fqnames=($(echo "${fq[@]%_*}" | tr ' ' '\n' | sort -u | tr '\n' ' '))

for i in ${fqnames[@]}; do

    echo "===================================="
    echo "Kallisto using ${SPECIES} cDNA for ${i}."
    ${KALLISTO_PATH}/kallisto quant \
                    -t ${CORENUM} \
                    -i ${REF_PATH}/ath.kindex \
                    -o ${ALIGN_PATH}/${i}_${SPECIES}_kallisto \
                    ${CLEAN_PATH}/${i}_R1.fq.gz ${CLEAN_PATH}/${i}_R2.fq.gz

    echo "HISAT2 using ${SPECIES} genome for ${i}."
    ${HISAT2_PATH}/hisat2 -p ${CORENUM} \
                  -x ${REF_PATH}/athht2index/genome \
                  -1 ${CLEAN_PATH}/${i}_R1.fq.gz \
                  -2 ${CLEAN_PATH}/${i}_R2.fq.gz \
                  -S ${ALIGN_PATH}/${i}_${SPECIES}_hisat2.sam

    ${SAMTOOL_PATH}/samtools sort \
                   -@ ${CORENUM} \
                   -o ${ALIGN_PATH}/${i}_${SPECIES}_hisat2.bam \
                   ${ALIGN_PATH}/${i}_${SPECIES}_hisat2.sam

    ${SAMTOOL_PATH}/samtools index ${ALIGN_PATH}/${i}_${SPECIES}_hisat2.bam

    ${RM_PATH} ${ALIGN_PATH}/${i}_${SPECIES}_hisat2.sam
    echo "====================================="

done
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~1stadd virus~~~~~~~~~~~~~~~~~~~~~~~~~~
fq=($(ls ${CLEAN_PATH} | grep .fq.gz))
fqnames=($(echo "${fq[@]%_*}" | tr ' ' '\n' | sort -u | tr '\n' ' '))

for i in ${fqnames[@]}; do

    echo "===================================="
    echo "Kallisto using ${SPECIES} cDNA for ${i}."
    ${KALLISTO_PATH}/kallisto quant \
                    -t ${CORENUM} \
                    -i ${REF_PATH}/ath_latentvirus1.kindex \
                    -o ${ALIGN_PATH}/${i}_${SPECIES}_kallisto \
                    ${CLEAN_PATH}/${i}_R1.fq.gz ${CLEAN_PATH}/${i}_R2.fq.gz

    echo "HISAT2 using ${SPECIES} genome for ${i}."
    ${HISAT2_PATH}/hisat2 -p ${CORENUM} \
                  -x ${REF_PATH}/athlatenvirus1ht2index/genome \
                  -1 ${CLEAN_PATH}/${i}_R1.fq.gz \
                  -2 ${CLEAN_PATH}/${i}_R2.fq.gz \
                  -S ${ALIGN_PATH}/${i}_${SPECIES}_hisat2.sam

    ${SAMTOOL_PATH}/samtools sort \
                   -@ ${CORENUM} \
                   -o ${ALIGN_PATH}/${i}_${SPECIES}_hisat2.bam \
                   ${ALIGN_PATH}/${i}_${SPECIES}_hisat2.sam

    ${SAMTOOL_PATH}/samtools index ${ALIGN_PATH}/${i}_${SPECIES}_hisat2.bam

    ${RM_PATH} ${ALIGN_PATH}/${i}_${SPECIES}_hisat2.sam
    echo "====================================="

done
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
