##
## originally by Yulong Niu
## yniu@mpipz.mpg.de
##

date

BIN_PATH=/usr/bin
SRATOOLS_PATH=/extDisk1/Biotools/sratoolkit.2.10.0-centos_linux64/bin
FASTP_PATH=/home/yniu/Biotools

RAW_PATH=/netscratch/dep_psl/grp_rgo/yniu/Flg22Acute/raw_data
CLEAN_PATH=/netscratch/dep_psl/grp_rgo/yniu/Flg22Acute/clean_data

CORENUM=40

#########################download SRA################################
cd ${RAW_PATH}

##~~~~~~~~~~~~~~~~~~~~single-end~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ref: https://www.nature.com/articles/nature21417
## repo: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA412447/
## ref: https://web.archive.org/web/20101114051536/http://codesnippets.joyent.com/posts/show/1826

SRARaw=($(seq 6983 1 7045))
SRARaw=("${SRARaw[@]/#/SRR623}")

SRAName=("flg22Pa_1h_1" "flg22417_1h_3" "flg22417_1h_2"
         "flg22417_1h_1" "WCS417_1h_3" "WCS417_1h_2"
         "WCS417_1h_1" "Mock_1h_3" "flg22Pa_1h_3"
         "flg22Pa_1h_2" "WCS417_0p5h_1" "WCS417_0p5h_2"
         "Mock_0p5h_2" "Mock_0h_3" "Mock_0p5h_1"
         "Mock_0h_1" "Mock_0p5h_3" "Mock_0h_2"
         "WCS417_0p5h_3" "flg22417_0p5h_1" "WCS417_6h_2"
         "WS417_6h_3" "Mock_6h_3" "WCS417_6h_1"
         "flg22417_6h_3" "flg22Pa_6h_1" "Mock_1h_1"
         "Mock_1h_2" "chitin_0p5h_2" "chitin_0p5h_3"
         "flg22Pa_0p5h_3" "chitin_0p5h_1" "flg22Pa_0p5h_1"
         "flg22Pa_0p5h_2" "flg22417_0p5h_2" "flg22417_0p5h_3"
         "chitin_1h_1" "chitin_1h_2" "chitin_1h_3"
         "Mock_3h_1" "Mock_3h_2" "Mock_3h_3"
         "WCS417_3h_1" "WCS417_3h_2" "WCS417_3h_3"
         "flg22417_3h_1" "flg22417_6h_2" "flg22417_6h_1"
         "chitin_6h_3" "flg22Pa_6h_2" "flg22Pa_6h_3"
         "flg22417_3h_3" "flg22417_3h_2" "flg22Pa_3h_2"
         "flg22Pa_3h_1" "chitin_3h_1" "flg22Pa_3h_3"
         "chitin_3h_3" "chitin_3h_2" "Mock_6h_2"
         "Mock_6h_1" "chitin_6h_2" "chitin_6h_1")

## check
IFS=$'\n' sorted=($(sort <<<"${SRAName[*]}")); unset IFS
printf "[%s]\n" "${sorted[@]}"

## for i in ${!SRARaw[*]}; do
mfile=(13 {29..31} 33 {35..38} 46 {54..56} 62)
for idx in ${!mfile[*]}; do

    i=${mfile[idx]}

    ${SRATOOLS_PATH}/fasterq-dump ${SRARaw[i]} -e ${CORENUM} \
                    -O ${RAW_PATH}

    ${BIN_PATH}/gzip ${SRARaw[i]}.fastq

    ${BIN_PATH}/mv ${SRARaw[i]}.fastq.gz ${SRAName[i]}.fq.gz

done
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#######################################################################

############################trim#######################################
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~raw~~~~~~~~~~~~~~~~~~
## sample names
fq=($(ls | grep .fq.gz))

for i in ${fq[@]}
do
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
#######################################################################
