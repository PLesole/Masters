# R script for data visualisation 

# load libraries
library(tidyverse)
library(ggvenn)
library(rtracklayer)
library(reshape2)
library(DEXSeq)
library(DRIMSeq)

# set working directory
setwd("C:/Palesa/School/Masters")

# import gtf file for data analysis
gtf <- import.gff("gencode.v44.primary_assembly.annotation.gtf")
gtf_df <- as.data.frame(gtf)
gtf_file <- unique(gtf_df[,c("gene_id","gene_name")])
saveRDS(gtf_file, "gtf.rds") # save for future use 

# load output files for each method
DEU <- readRDS('DEU_results.rds')
DTU <- readRDS('DTU_dmDS.rds')
sigDTU <- readRDS('DTU_stageR_results.rds')
DASE <- readRDS('DASE_results.rds')

#### define function for plotting ####
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

#### Alignment rates ####
genomic_rates <- data.frame(Sample = c("DMSO_rep1","DMSO_rep2","DMSO_rep3","PMA_rep1","PMA_rep2","PMA_rep3"),
                            Rate = as.numeric(c("89.24","87.30","88.01","85.67","88.37","87.85")))
genomic_rates<- mutate(genomic_rates, Condition = case_when(grepl("PMA", Sample) ~ "Macrophages",grepl("DMSO", Sample) ~ "Monocytes"))


transcriptomic_rates <- data.frame(Sample = c("DMSO_rep1","DMSO_rep2","DMSO_rep3","PMA_rep1","PMA_rep2","PMA_rep3"),
                                   Rate = as.numeric(c("82.16","79.75","80.87","78.82","81.54","81.34")))
transcriptomic_rates <- mutate(transcriptomic_rates, Condition = case_when(grepl("PMA", Sample) ~ "Macrophages",grepl("DMSO", Sample) ~ "Monocytes"))


genome_rate <- ggplot(genomic_rates, aes(x = Sample, y = Rate, fill = Condition)) + geom_bar(stat = "identity") +
  coord_flip() +
  scale_y_continuous(expand = c(0,0)) + 
  scale_fill_manual(values = c("#46B1E1","#78206E")) +
  labs(x = "Sample", y = "Genomic alignment rate (%)", fill = "Condition") +
  plot_theme() # plot genomic rates 


transcriptome_rate <- ggplot(transcriptomic_rates, aes(x = Sample, y = Rate, fill = Condition)) + geom_bar(stat = "identity") +
  coord_flip() +
  scale_y_continuous(expand = c(0,0)) + 
  scale_fill_manual(values = c("#46B1E1","#78206E")) +
  labs(x = "Sample", y = "Transcriptomic alignment rate (%)", fill = "Condition") +
  plot_theme() # plot transcriptomic rates 

ggsave("Genomic_alignment_rates.png", genome_rate, units = "mm", height = 80, width = 120, dpi = 500)
ggsave("Transcriptomic_alignment_rates.png", transcriptome_rate, units = "mm", height = 80, width = 120, dpi = 500)

#### FastQC - per base sequence quality ####
options(scipen = 999) #collapse exponential form 

per_base_quality <- read_tsv("multiqc_dmso/fastqc_per_base_sequence_quality_plot.tsv")

colnames(per_base_quality) <- c("Position (bp)", paste0(rep(c("DMSO","PMA"), each = 6), "_rep", rep(1:3, each = 2), "_R", rep(1:2)))

per_base_quality_long <- pivot_longer(per_base_quality, cols = -c("Position (bp)"), 
                                          names_to = 'Sample', values_to = 'Quality_scores')

per_base_quality_long <- mutate(per_base_quality_long, Condition = case_when(grepl("DMSO", Sample) ~ "Monocytes",
                                                                                     grepl("PMA", Sample)  ~ "Macrophages"))

per_base_quality_long <- mutate(per_base_quality_long, Read_num = case_when(grepl("R1",Sample) ~ "R1", grepl("R2",Sample) ~ "R2")) # assign read numbers 
per_base_quality_long$Read_num <- as.factor(per_base_quality_long$Read_num)
per_base_quality_long$Condition <- factor(per_base_quality_long$Condition, levels = c('Monocytes', 'Macrophages'))

per_base_plot <- ggplot(per_base_quality_pma_long, aes(x = `Position (bp)`, y = Quality_scores, color = Condition)) + geom_line(aes(group = Sample), linewidth = 0.7) + 
  facet_wrap(~Read_num) +
  labs(x = "Base Position (bp)", y = "Quality score") +
  scale_colour_manual(values = c("Monocytes" = "#46B1E1",
                                 "Macrophages" = "#78206E")) + 
  plot_theme()

ggsave("per_base_sequence_quality.png", per_base_plot, units = "mm", height = 90 , width = 120 , dpi = 500)

#### PCA ####
##### Exon level #### 
# extract normalised counts from results dataSet
DEU_counts <- data.frame(DEXSeq::counts(DEU, normalized = TRUE))
DEU_counts <- t(DEU_counts) # flip rows and columns 

zero_var <- which(apply(DEU_counts, 2, var) == 0) # remove columns with zero variance
DEU_counts_filtered <- DEU_counts[, -zero_var]

pca_res <- prcomp(DEU_counts_filtered, scale. = TRUE) # calculate PCs
percentVar <- pca_res$sdev^2 / sum(pca_res$sdev^2) * 100 # Percentage of variance explained by each component

pca_df <- data.frame(PC1 = pca_res$x[, 1], 
                     PC2 = pca_res$x[, 2],
                     Condition = factor(rep(c("Monocytes", "Macrophages"), each = 3))) # Create a data frame for ggplot

# Create the PCA plot using ggplot2
exon_PCA <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", round(percentVar[1], 2), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2], 2), "% variance")) +
  scale_color_manual(values = c("Monocytes" = "#46B1E1", "Macrophages" = "#78206E")) + 
  plot_theme()

ggsave("exon_level_PCA.png", exon_PCA, units = "mm", height = 80, width = 100, dpi = 500)

##### Transcript level ####
# extract transcript ratios from dmDS
DTU_counts <- DRIMSeq::counts(DTU)
rownames(DTU_counts) <- paste0(DTU_counts[,"gene_id"], ":", DTU_counts[,"feature_id"]) # create row names of gene and transcript IDs
DTU_counts <- DTU_counts[, -c(1,2)] # remove extra columns 
DTU_counts <- t(DTU_counts) # flip rows and columns 

pca_res_dtu <- prcomp(DTU_counts, scale. = TRUE) # calculate PCs
percentVar_dtu <- pca_res_dtu$sdev^2 / sum(pca_res_dtu$sdev^2) * 100 # Percentage of variance explained by each component

pca_dtu <- data.frame(PC1 = pca_res_dtu$x[, 1], 
                      PC2 = pca_res_dtu$x[, 2],
                      Condition = factor(rep(c("Monocytes", "Macrophages"), each = 3))) # Create a data frame for ggplot

# Create the PCA plot using ggplot2
transcript_PCA <- ggplot(pca_dtu, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", round(percentVar_dtu[1], 2), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar_dtu[2], 2), "% variance")) +
  scale_color_manual(values = c("Monocytes" = "#46B1E1", "Macrophages" = "#78206E")) + 
  plot_theme()

ggsave("transcript_level_PCA.png",transcript_PCA, units = "mm", height = 80, width = 100, dpi = 500)

##### Gene level aggregation - transcript ####
rownames(DTU_counts) <- DTU_counts[,"feature_id"] # create row names of gene and transcript IDs

tx2gene <- data.frame(transcript_id = rownames(DTU_counts),
                      gene_id = DTU_counts$gene_id) # map genes to transcripts to sum to genes

Gene_counts_dtu <- DTU_counts %>% mutate(gene_id = tx2gene$gene_id) %>%
  group_by(gene_id) %>% dplyr::summarise(across(where(is.numeric), sum))

rownames(Gene_counts_dtu) <- Gene_counts_dtu$gene_id # set genes as rownamess
Gene_counts_dtu <- Gene_counts_dtu[,-1] # remove gene IDs column

Gene_counts_dtu <- t(Gene_counts_dtu) # flip rows and columns 

pca_res_gt <- prcomp(Gene_counts_dtu, scale. = TRUE) # calculate PCs
percentVar_gt <- pca_res_gt$sdev^2 / sum(pca_res_gt$sdev^2) * 100 # Percentage of variance explained by each component

pca_gt <- data.frame(PC1 = pca_res_gt$x[, 1], 
                     PC2 = pca_res_gt$x[, 2],
                     Condition = factor(rep(c("Monocytes", "Macrophages"), each = 3))) # Create a data frame for ggplot

# Create the PCA plot using ggplot2
genetrans_PCA <- ggplot(pca_gt, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", round(percentVar_gt[1], 2), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar_gt[2], 2), "% variance")) +
  scale_color_manual(values = c("Monocytes" = "#46B1E1", "Macrophages" = "#78206E")) + 
  plot_theme()

ggsave("gene_aggregate_transcript_PCA.png", genetrans_PCA, units = "mm", height = 80, width = 100, dpi = 500)

##### Gene level aggregation - exon ####
DEU_counts$gene_id <- sapply(strsplit(rownames(DEU_counts), ":"), `[`,1) # create gene id using rownames
Gene_counts_deu <- DEU_counts %>% group_by(gene_id) %>% summarise(across(where(is.numeric),sum))
rownames(Gene_counts_deu) <- Gene_counts_deu$gene_id # set gene ids as rownames 

Gene_counts_deu <- Gene_counts_deu[,-1] # remove gene IDs
Gene_counts_deu <- t(Gene_counts_deu) # flip rows and columns 

zero_var_ge <- which(apply(Gene_counts_deu, 2, var) == 0) # remove columns with zero variance
Gene_counts_deu_filtered <- Gene_counts_deu[, -zero_var_ge]

pca_res_ge <- prcomp(Gene_counts_deu_filtered, scale. = TRUE) # calculate PCs
percentVar_ge <- pca_res_ge$sdev^2 / sum(pca_res_ge$sdev^2) * 100 # Percentage of variance explained by each component

pca_ge <- data.frame(PC1 = pca_res_ge$x[, 1], 
                     PC2 = pca_res_ge$x[, 2],
                     Condition = factor(rep(c("Monocytes", "Macrophages"), each = 3))) # Create a data frame for ggplot

# Create the PCA plot using ggplot2
geneexons_PCA <- ggplot(pca_ge, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", round(percentVar_ge[1], 2), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar_ge[2], 2), "% variance")) +
  scale_color_manual(values = c("Monocytes" = "#46B1E1", "Macrophages" = "#78206E")) + 
  plot_theme()

ggsave("gene_aggregate_exon_PCA.png", geneexons_PCA, units = "mm", height = 80, width = 100, dpi = 500)

#### Number of differential splicing events and affected genes ####
# Create a data frame with feature counts and gene count for each differentiation protocol
event_counts <- data.frame(
  Method = c("DEU", "DTU", "DAS"),
  Significant_Features = c(nrow(DEU), nrow(sigDTU), nrow(DASE)),
  Affected_Genes = c(length(unique(DEU$gene_id)),length(unique(sigDTU$geneID)),length(unique(DASE$gene_id))
  )) 

events_df <- pivot_longer(event_counts, cols = c("Significant_Features","Affected_Genes"),
                          names_to = "Category", values_to = "Counts")

# filter by count-based method
deu_events <- filter(events_df, Method == "DEU")
dtu_events <- filter(events_df, Method == "DTU")
das_events <- filter(events_df, Method == "DAS")


deu_plot <- ggplot(deu_events, aes(x = Method, y = Counts, fill = Category)) + geom_bar(stat = "identity", position = "dodge") +
  coord_flip() + 
  scale_fill_manual(values = c("#46B1E1","#78206E"), labels = c("Affected genes","Number detected")) + 
  scale_y_continuous(expand = c(0,0)) +  
  labs(x = "Method", y = "Count", fill = "Category") + 
  plot_theme()

dtu_plot <- ggplot(dtu_events, aes(x = Method, y = Counts, fill = Category)) + geom_bar(stat = "identity", position = "dodge") +
  coord_flip() + 
  scale_fill_manual(values = c("#46B1E1","#78206E"), labels = c("Affected genes","Number detected")) + 
  scale_y_continuous(expand = c(0,0)) +  
  labs(x = "Method", y = "Count", fill = "Category") + 
  plot_theme()

das_plot <- ggplot(das_events, aes(x = Method, y = Counts, fill = Category)) + geom_bar(stat = "identity", position = "dodge") +
  coord_flip() + 
  scale_fill_manual(values = c("#46B1E1","#78206E"), labels = c("Affected genes","Number detected")) + 
  scale_y_continuous(expand = c(0,0)) +  
  labs(x = "Method", y = "Count", fill = "Category") + 
  plot_theme()

ggsave("DEU_Events.png", deu_plot, units = "mm", height = 50, width = 100, dpi = 500)
ggsave("DTU_Events.png", dtu_plot, units = "mm", height = 50, width = 100, dpi = 500)
ggsave("DAS_Events.png", das_plot, units = "mm", height = 50, width = 100, dpi = 500)


#### Overlap between differentially spliced genes ####
# pull  out unique differentially spliced geneIDs
DEU_genes <- unique(DEU$gene_id)
DTU_genes <- unique(sigDTU$geneID)
DAS_genes <- unique(DAS$gene_id)

# pull out gene names
DEU_name <- gtf_file[gtf_file$gene_id %in% DEU_genes,]
DTU_name <- gtf_file[gtf_file$gene_id %in% DTU_genes,]
DAS_name <- gtf_file[gtf_file$gene_id %in% DAS_genes,]

genes_overlap <- ggvenn(list(DEU = DEU_genes, DTU = DTU_genes, DAS = DAS_genes), 
                      fill_color = c("#46B1E1", "#F4A259","#78206E"),
                      stroke_color = "lightgrey",
                      stroke_size = 0.01, 
                      show_percentage = FALSE,
                      text_size = 7) 

ggsave("genes_Overlap.png", genes_overlap, units = "mm", dpi = 500)
