#!/bin/bash
# Shell script to analysis paired-end RNA-seq reads using FASTQC and multiqc

Input_dir=~/Palesa/Masters/Green_data
Output_dir=~/Palesa/Masters/Green_data/fastqc

#check if inout and output directories exist; make output directory if it does not exist
if [[ -z $Input_folder ]]; then
	echo -e "The input directory is missing"
	exit -1
fi

if [[ -z $Output_dir ]]; then
	echo -e "The output directory is missing"
  echo -e "Making output directory now"
  mkdir -p "$Output_dir"
	exit -1
fi



for file in $(ls ${Input_dir}/*.fastq.gz  | sed 's/_[12].fastq.gz$//' | uniq)
  do

    fastqc ${file}_1.fastq.gz ${file}_2.fastqc.gz --outdir $Output_dir
  done

multiqc $Output_dir
