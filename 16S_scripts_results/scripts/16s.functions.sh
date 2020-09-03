#!/bin/bash

# scripts for 16S data analysis
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

log() {

    echo $(date -u)": "$1 >> $logfile

}

cleanup() {
 
    find $working_dir -type f -print0 | xargs -0 rm -f

}

sampleSizes() {
    
    local seqs=$1
    local mapping_file=$2

    for sample in $(cut -f 1 $mapping_file | tail -n +2)
    do 

        paste <(echo $sample) <(grep -c ">"$sample"_" $seqs)

    done | sort -n -k2,2

}

