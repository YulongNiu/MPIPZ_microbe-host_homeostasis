#!/bin/bash

# scripts for 16S data analysis
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

# exits whenever a function returns 1
set -e

# get path to scripts
scripts_dir=$(dirname $0)

# load config file
config_file=$1
source $config_file

# load functions
source $scripts_dir/16s.functions.sh

# activate QIIME, etc.
source $scripts_dir/activate.sh

# cleanup
### cleanup

# processing of the raw reads for all libraries:
# joining of pair reads, quality control, demultiplexing

if [ "$l_list" != "" ]
then

    rm -f $working_dir/seqs.fasta

    for l in $l_list
    do

       rm -f -r $working_dir/"$l"
       mkdir -p $working_dir/"$l"

       # merge paired reads
       log "["$l"] joining paired-end reads..."
       add_barcode_to_label.py $data_dir/"$l"_forward_reads.fastq.gz \
                               $data_dir/"$l"_barcodes.fastq.gz \
                               $working_dir/"$l"/forward_reads.fastq \
                               &>> $output
       add_barcode_to_label.py $data_dir/"$l"_reverse_reads.fastq.gz \
                               $data_dir/"$l"_barcodes.fastq.gz \
                               $working_dir/"$l"/reverse_reads.fastq \
                               &>> $output
       usearch -fastq_mergepairs $working_dir/"$l"/forward_reads.fastq \
               -reverse $working_dir/"$l"/reverse_reads.fastq \
               -fastqout $working_dir/"$l"/joined.fastq \
               &>> $output
       extract_barcodes.py -f $working_dir/"$l"/joined.fastq \
                           -c barcode_in_label \
                           --char_delineator 'BC=' \
                           --bc1_len $bc_length \
                           -o $working_dir/"$l"/ \
                           &>> $output

        # demultiplex
        log "["$l"] demultiplexing and quality filtering..."
        split_libraries_fastq.py -i $working_dir/$l/joined.fastq \
                                 -b $working_dir/$l/barcodes.fastq \
                                 --rev_comp_mapping_barcodes \
                                 --max_barcode_errors $bc_err \
                                 --barcode_type $bc_length \
                                 -o $working_dir/$l \
                                 -m $data_dir/"$l"_mapping.txt \
                                 -q $phred \
                                 --phred_offset 33 \
                                 &>> $output

        # save and edit barcode label identifier for usearch compatibility
        rm -f $working_dir/$l/seqs.fasta
        cat $working_dir/$l/seqs.fna | \
            awk '/>/ {print $0, "barcodelabel="$1} !/>/ {print $0}' | \
            sed 's/=>/=/g;s/_[0-9]*$/;/g;s/ /;/g' >> \
            $working_dir/$l/seqs.fasta

        # concatenate all demultiplexed sequences
        cat $working_dir/$l/seqs.fasta >> $working_dir/seqs.fasta

        mv $working_dir/"$l"/split_library_log.txt $working_dir/"$l"_split_library_log.txt

        ### rm -f -r $working_dir/"$l"

    done

fi

# reference-based OTU clustering using the UPARSE-REF algorithm
log "reference-based OTU clustering using UPARSE-REF..."
usearch -uparse_ref $working_dir/seqs.fasta \
        -db $ref_seqs \
        -strand plus \
        -uparseout $working_dir/otu_clustering.up \
        -threads $n_cores \
        &>> $output

# generate UC-compatible mapping file for USEARCH compatibility
rm -f $working_dir/otu_clustering.uc
cat $working_dir/otu_clustering.up | \
    awk -v i="$id_threshold" '{if ($3<i*100 || $2=="chimera" || $2=="other") {$5=$2}
                               print $1, $5}' | \
    sed 's/ /\t/g;s/^/H\t\t\t\t\t\t\t\t/g' \
    1>> $working_dir/otu_clustering.uc \
    2>> $output

# convert UC file to txt
log "converting UC OTU table file into text format..."
python $usearch_dir/uc2otutab.py $working_dir/otu_clustering.uc \
    1> $working_dir/otu_table_unfiltered.txt \
    2>> $output

# remove chimera and other entries from OTU table
# NOTE: the unfiltered OTU table is useful to check
# the percentage of posible contamination ('other')
# or of chimeric sequences ('chimera') per sample
grep -v -e other \
     -v -e chimera \
     -v -e noisy \
     -v -e good \
     $working_dir/otu_table_unfiltered.txt \
     1> $working_dir/otu_table.txt \
     2> $output

# TODO:
# total sum normalization
# calculate alpha and beta diversity indices

# parse OTU TXT tables

sed -i 's/OTUId.//g' $working_dir/otu_table.txt
sed -i 's/OTUId.//g' $working_dir/otu_table_unfiltered.txt

log "DONE!"
