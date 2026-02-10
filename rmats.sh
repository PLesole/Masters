#!/bin/bash

wdir=~/Palesa/Masters 
GTF=gencode.v44.primary_assembly.annotation.gtf
read_length=150

#code to run prep and post step simultaneously
 python rmats.py --b1 CNT_bam.txt --b2 PMA_bam.txt 
    --gtf $wdir/$GTF -t paired --readLength $read_length --paired-stats
    --nthread 8 --task both --od $wdir --tmp $wdir 
# time print out time taken to run code
# --b1 and b2 text file containing paths to bam files for sample 1 and 2
# --gtf reference gtf used in read alignment 
# -t indicates library type
# --task both runs both prep and prep step 
# output directory (od) and directory for tmp files (tmp) in current working directory     
# required read length to determine intron length for PSI calculations 
# --paired-stats toggles on paired statistical model that runs within PAIRADAISE
# nthreads indicates number of threads used to run -- task          

# code to run stat step
 python rmats.py --od $wdir --tmp $wdir --task stat 
# od and tmp must be directory containing tmp and output files after post step  