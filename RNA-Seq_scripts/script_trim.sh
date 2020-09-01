
## originally by Yulong Niu
## yulong.niu@hotmail.com

date

RAW_PATH=/netscratch/dep_psl/grp_rgo/yniu/KaWaiFlg22/raw_data_soil
CLEAN_PATH=/netscratch/dep_psl/grp_rgo/yniu/KaWaiFlg22/clean_data_soil

FASTP_PATH=/home/yniu/Biotools
FASTQC_PATH=/opt/share/software/bin

CORENUM=16

cd ${RAW_PATH}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~raw~~~~~~~~~~~~~~~~~~
## sample names
fq=($(ls | grep .fq.gz))

for i in ${fq[@]}
do
    echo "FastQC ${i}"
    ${FASTQC_PATH}/fastqc -o ${RAW_PATH} \
                  -t ${CORENUM} \
                  ${i}

    echo "Trimming ${i}."
    ${FASTP_PATH}/fastp -w ${CORENUM} \
                 -z 6 \
                 -p \
                 -U --umi_loc=read1 --umi_len=8 \
                 --trim_tail1=2 \
                 -h ${i%%.*}.html \
                 -i ${i} \
                 -o ${CLEAN_PATH}/${i}
done
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~1stadd~~~~~~~~~~~~~~~~~~~~~~~~~~
fq=($(ls | grep fq.gz))
fqnames=($(echo "${fq[@]%_*}" | tr ' ' '\n' | sort -u | tr '\n' ' '))

for i in ${fqnames[@]}
do
    echo "Trimming ${i}_R1.fq.gz ${i}_R2.fq.gz."

    ${FASTP_PATH}/fastp -w ${CORENUM} \
                 -z 6 \
                 -p -c \
                 -h ${i}.html \
                 -i ${i}_R1.fq.gz -I ${i}_R2.fq.gz \
                 -o ${CLEAN_PATH}/${i}_R1.fq.gz -O ${CLEAN_PATH}/${i}_R2.fq.gz
done
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~soil~~~~~~~~~~~~~~~~~~~~~~~~~~
fq=($(ls | grep fq.gz))
fqnames=($(echo "${fq[@]%_*}" | tr ' ' '\n' | sort -u | tr '\n' ' '))

for i in ${fqnames[@]}
do
    echo "Trimming ${i}_R1.fq.gz ${i}_R2.fq.gz."

    ${FASTP_PATH}/fastp -w ${CORENUM} \
                 -z 6 \
                 -p -c \
                 -h ${i}.html \
                 -i ${i}_R1.fq.gz -I ${i}_R2.fq.gz \
                 -o ${CLEAN_PATH}/${i}_R1.fq.gz -O ${CLEAN_PATH}/${i}_R2.fq.gz
done
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

date
