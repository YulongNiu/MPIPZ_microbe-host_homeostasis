date

RAW_PATH=/netscratch/dep_psl/grp_rgo/yniu/KaWaiFlg22/raw_data
CLEAN_PATH=/netscratch/dep_psl/grp_rgo/yniu/KaWaiFlg22/clean_data

FASTP_PATH=/home/yniu/Biotools
FASTQC_PATH=/opt/share/software/bin

CORENUM=16

cd ${RAW_PATH}

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

date
