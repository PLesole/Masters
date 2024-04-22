#!/bin/bash

Input_folder=~/Palesa/Masters/Green_data
Output_folder=~/Palesa/Masters/Green_data/fastqc

mkdir -p "$Output_folder"


for file in SRR8926873 SRR8926874 SRR8926875
  do
    fastq_R1="${Input_folder}/${file}_1.fastq.gz"
    fastq_R2="${Input_folder}/${file}_2.fastq.gz"
    fastqc $fastq_R1 $fastq_R2 --outdir $Output_folder
  done

multiqc $Output_folder
