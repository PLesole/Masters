#!/bin/bash

# bash script to perform transcriptome alignment using RSEM and HISAT2

#### Prepare reference sequences that will be used to perform transcriptome alignment in HISAT2 ####
Genome_fasta=~/Palesa/Masters/GRCh38.primary_assembly.fa
GTF=~/Palesa/Masters/gencode.v44.primary_assembly.annotation.gtf
HISAT2=~/Palesa/hisat2-2.2.1-Linux_x86_64/hisat2-2.2.1
rsem=~/Palesa/Masters/RSEM/RSEM-1.3.3
ref_name=GRCh38
threads=6

$rsem/rsem-prepare-reference --gtf $GTF -p $threads --hisat2-hca --hisat2-path $HISAT2 $Genome_fasta $ref_name

### calculate transcript abundance ####
wdir=~/Palesa/Masters/RSEM
input_dir=~/Palesa/Masters/Data

for file in $(ls ${input_dir}/*.fastq.gz  | sed 's/_[12].fastq.gz$//' | uniq)
 do
       # Testing if the two paired files exist
       if [[ -z ${file}_1.fastq.gz  ]] || [[ -z ${file}_2.fastq.gz ]]; then
           echo -e "Read 1 and Read 2 file cannot be detected for $file"
           exit -1
       fi
       output=$(basename $file)
    $rsem/rsem-calculate-expression --hisat2-hca --hisat2-path $HISAT2 --paired-end -p $threads --strandedness reverse --no-bam-output ${file}_R1_001.fastq.gz ${file}_R2_001.fastq.gz $ref_name $wdir/$output
 done
