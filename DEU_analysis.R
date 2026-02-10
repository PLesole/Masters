## R script for DEU analysis using DEXSeq and quality control metrics (exome complexity analysis and exon read count distribution)

# load libraries 
library(DEXSeq)
library(tidyverse)
library(GenomicFeatures)
library(txdbmaker) 
library(Rsamtools)
library(BiocParallel)
library(GenomicAlignments)
library(ggpubr)

# set working directory 
setwd("~/Palesa/Masters")

### define function for plotting ####
plot_theme <- function() {
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "grey20"),
        axis.text.x = element_text(
          size = 11, 
          colour = "black"
        ),
        axis.text.y = element_text(
          size = 11, 
          colour = "black"
        ),
        axis.title.x = element_text(
          size = 13, 
          colour = "black"
        ),
        axis.title.y = element_text(
          size = 13, 
          colour = "black"
        ))
}

#### DEU analysis ####

## create flattened annotation ##
txdb <- makeTxDbFromGFF("~/Palesa/Masters/gencode.v44.primary_assembly.annotation.gtf") #create TxDb from gtf file used to get alignments and get flattened annotation 
flat <- exonicParts(txdb, linked.to.single.gene.only = TRUE) # extract disjointed exonic parts from annotation and ensures that they are no not overlap into other genes
names(flat) <- sprintf("%s:E%0.3d", flat$gene_id, flat$exonic_part) # annotate exonic parts

# create object containing paths to BAM files
bamfiles = c("~/Palesa/Masters/BAM_files/CNT_replicate_1.bam", 
             "~/Palesa/Masters/BAM_files/CNT_replicate_2.bam",
             "~/Palesa/Masters/BAM_files/CNT_replicate_3.bam",
             "~/Palesa/Masters/BAM_files/PMA_replicate_1_sorted.bam",
             "~/Palesa/Masters/BAM_files/PMA_replicate_2_sorted.bam",
             "~/Palesa/Masters/BAM_files/PMA_replicate_3_sorted.bam")
bamfiles <- BamFileList(bamfiles)

param <- MulticoreParam(workers = multicoreWorkers()) # parallel to work on multiple cores

# count reads using summarizedOverlaps and create summarizedExperiment from flattened annotation and read counts 
se <- summarizeOverlaps(flat, BamFileList(bamfiles), singleEnd = FALSE,
                        fragments = TRUE, ignore.strand = FALSE, BPPARAM = param) # ensure that counts no not ignore library type and strandedness

colData(se)$condition <- factor(c("CNT", "CNT", "CNT", "PMA", "PMA", "PMA")) # specify treatment conditions: CNT - Control, PMA - Treated 
colData(se)$libType <- "Paired-end" # specify library type - reads are from a paired-end experiment 

dxd <- DEXSeqDataSetFromSE(se, design = ~sample + exon + condition:exon) # create DEXSeqDataset from summarizedExperiment

dxd <- estimateSizeFactors(dxd) # data normalisation according to DESeq
dxd <- estimateDispersions(dxd) # estimate dispersions (variability) of data
dxd <- testForDEU(dxd) # test for exon usage between conditions by fitting model design on genes 
dxd <- estimateExonFoldChanges(dxd) # estimate relative exon usage fold change using coefficients of model

dxr <- DEXSeqResults(dxd, independentFiltering = TRUE) # extract results of analysis with independent filtering of significant exons 

saveRDS(dxr, 'DEU_results.rds') # save RDS for downstream analysis

### exon read count distribution ####
exon_counts <- counts(dxr, normalized = TRUE) # extract normalised exon counts

# filter out lowly expressed exons
keep_exons <- rowMeans(exon_counts >=5) 
exon_counts_filter <- exon_counts[keep_exons, ]

exon_long <- melt(log2(exon_counts_filter + 1), varnames=c("exon","sample"), value.name="logCounts")

exon_long$condition <- factor(ifelse(grepl("CNT", exon_long$sample), "Monocyte", "Macrophage"),
                              levels = c("Monocyte","Macrophage")) # change factors to cell type 


deu_counts_violin <- ggplot(exon_long, aes(x = condition, y = logCounts, fill = condition)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
  stat_summary(fun = "median", geom = "point", size = 2, color = "black") +
  labs(x = "Cell Type", y = "log2(Counts)") +
  scale_fill_manual(values = c("Monocyte" = "#46B1E1", "Macrophage" = "#78206E")) +
  plot_theme() ## making violin plots to show expression bias 

ggsave("exon_counts_distribution.png", deu_counts_violin, units = "mm", height = 80, width = 100 , dpi = 500)


#### exome complexity analysis ####
exon_df <- as.data.frame(rowRanges(se)) # convert summarizedExperiment into data frame 

# Summarize counts per gene and remove NA values
exon_complexity <- as.data.frame(assay(se)) %>%
  mutate(gene_id = exon_df$gene_id) %>%
  group_by(gene_id) %>%
  summarise(n_exons = n(),
            mean_expr = mean(rowMeans(across(where(is.numeric)))),
            var_expr = var(rowMeans(across(where(is.numeric))))) %>%
  na.omit()

exome_complexity <- ggplot(exon_complexity, aes(x = n_exons, y = var_expr)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  scale_x_continuous(breaks = pretty(exon_complexity$n_exons)) +
  labs(x = "Number of exons per gene", y = "Variance in exon-level expression",
       title = "Exon usage complexity per gene (DEXSeq input)") +
  plot_theme()

ggsave("exome_complexity.png", exome_complexity, units = "mm" , dpi = 500)