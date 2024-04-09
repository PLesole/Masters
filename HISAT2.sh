#!/bin/bash

REF_GENOME_INDEX="~/Palesa/Masters/hisat2_v44_genome_tran/hisat2_v44_genome_tran"
OD=~/Palesa/Masters

for file in UNT_replicate_1 UNT_replicate_2 UNT_replicate_3
 do
    replicate_R1=${file}_R1_001.fastq.gz
    replicate_R2=${file}_R2_001.fastq.gz
    output_prefix=$file
    ~/Palesa/hisat2-2.2.1-Linux_x86_64/hisat2-2.2.1/hisat2 -x $REF_GENOME_INDEX -1 $replicate_R1 -2 $replicate_R2 --no-softclip -S $OD/${output_prefix}.sam
 done

for file in PMA_replicate_1 PMA_replicate_2 PMA_replicate_3
do
    replicate_R1=${file}_R1_001.fastq.gz
    replicate_R2=${file}_R2_001.fastq.gz
    output_prefix=$file
    ~/Palesa/hisat2-2.2.1-Linux_x86_64/hisat2-2.2.1/hisat2 -x $REF_GENOME_INDEX -1 $replicate_R1 -2 $replicate_R2 --no-softclip -S $OD/${output_prefix}.sam
done
