#!/bin/bash

# scripts for 16S data analysis
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

# list of library IDs, separated by single spaces
l_list="werewolf_001"

# working_dir is the folder where (intermediate) results are stored
working_dir="/netscratch/dep_psl/grp_rgo/garridoo/werewolf_001/"

# data_dir is the folder where the raw data (and maping file) are deposited
data_dir="/biodata/dep_psl/grp_rgo/garridoo/werewolf/data/"
ref_seqs=$data_dir"/werewolf_001_ref_seqs.fasta"
mapping_file=$data_dir"/werewolf_001_mapping.txt"

# paths to relevant reference database files
refdata_dir="/biodata/dep_psl/grp_rgo/tools/16s/ref_data"
gg_core_db=$refdata_dir"/gg_13_8_otus/rep_set/97_otus.fasta"
gg_core_aligned_db=$refdata_dir"/gg_13_8_otus/rep_set_aligned/97_otus.fasta"
gg_taxonomy=$refdata_dir"/gg_13_8_otus/taxonomy/97_otu_taxonomy.txt"
gold_db=$refdata_dir"/cs_gold.fa"

# path to folder containing usearch binary and python scripts
usearch_dir="/biodata/dep_psl/grp_rgo/tools/usearch/"

### debug

output=$working_dir"/output.txt"
logfile=$working_dir"/log.txt"

### parameters

# number of threads for parallel steps
n_cores=30
# min Phred score for quality control
phred=30
# barcode length
bc_length=24
# max. barcode errors
bc_err=2
# min. number of reads for a given OTU
min_size=2
# identity threshold for OTU clustering
id_threshold=1
# min length of alingment required to sequences in the database
min_id_aln=0.75
# subsampling depth used for alpha-diversity
subsampling_depth=1500
# OTU clustering method
cls_method="uclust"

