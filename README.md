# Masters - Evaluating Count-Based Methods for Alternative Splicing Analysis in Monocyte-to-macrophage Differentiation 

This repository contains code for an end-to-end analysis pipeline to evaluate count-based alternative splicing analysis methods
## This repo contains:
  1. A shell script to perform QC on sequencing reads using FastQC and MultiQC
  2. A shell script to perform read alignment to the human reference genome using HISAT2
  3. A shell script to perform read alignment to the human transcriptome and transcript-level quantification using RSEM according to the Human Cell Atlas (HCA) Smart-Seq pipeline
  4. A R script to perform DEU analysis using DEXSeq and quality control - exome complexity analysis and exon read count distribution
  5. A R script to perform DTU analysis using DRIMSeq and stage-wise analysis using stageR and quality control - transcriptome complexity analysis and estimated transcript ratio distribution
  6. A shell script to perform DAS event analysis using rMATS-turbo
  7. A R script to perform DAS event analysis quality control - mean junction count distribution and filtering output from rMATS-turbo's prep step
  8. A R script for data visualisation - Alignment rates, PCA, per sequence quality (FastQC), Number of differentially splicing genes and overlap between differentially spliced genes
  9. A R script to create visualisation plots for representative genes uniquely identified by each count-based method (DEU and DTU analysis)
  10. A shell script to create sashimi plots using ggsashaimi within a Docker container
