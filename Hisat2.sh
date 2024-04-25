#!/bin/bash
# Shell script to align reads to reference genome

REF_GENOME_INDEX=~/Palesa/Masters/hisat2_v44_genome_tran/hisat2_v44_genome_tran
OD=~/Palesa/Masters/Green_data
WDIR=~/Palesa/hisat2-2.2.1-Linux_x86_64/hisat2-2.2.1

# make directories if not exist and enter working directory.
[[ ! -d ${OD} ]] && mkdir -p ${OD}
cd ${OD}

#check necessary parameters
if [[ -z $REF_GENOME_INDEX ]]; then
	echo -e "The index file for HISAT2 is missing"
	exit -1
fi

if [[ -z $OD ]]; then
	echo -e "The input directory for HISAT2 is missing"
	exit -1
fi

for file in $(ls ${OD}/*.fastq.gz  | sed 's/_[12].fastq.gz$//' | uniq)
 do
       # Testing if the two paired files exist
       if [[ -z ${file}_1.fastq.gz  ]] || [[ -z ${file}_2.fastq.gz ]]; then
           echo -e "Read 1 and Read 2 file cannot be detected for $file"
           exit -1
       fi
    output_prefix=$(basename $file)
    $WDIR/hisat2 -x $REF_GENOME_INDEX -1 ${file}_1.fastq.gz -2 ${file}_2.fastq.gz --no-softclip -S $OD/${output_prefix}.sam
 done
