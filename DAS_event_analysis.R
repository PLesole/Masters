# R script for quality control analysis from DAS event analysis using rMATS-turbo

# load libraries
library(tidyverse)
library(GenomicFeatures)
library(txdbmaker)

# set working directory
setwd("~/Palesa/Masters")

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

#### junction count distribution ####
# load output files from rMATS
SE <- read_tsv("DS_Analysis/CNTvsPMA/SE.MATS.JC.txt")
RI <- read_tsv("DS_Analysis/CNTvsPMA/RI.MATS.JC.txt")
A5SS <- read_tsv("DS_Analysis/CNTvsPMA/A3SS.MATS.JC.txt")
MXE <- read_tsv("DS_Analysis/CNTvsPMA/MXE.MATS.JC.txt")

# assign events to file type 
SE$ID...12 <- "SE"
MXE$ID...12 <- "MXE"
RI$ID...12 <- "RI"
A3SS$ID...12 <- "A3SS"
A5SS$ID...12 <- "A5SS"

# filter unimportant columns
SE_df <- dplyr::select(SE, ID...1, GeneID, geneSymbol, chr, ID...12,FDR, IncLevel1, IncLevel2, IncLevelDifference)
colnames(SE_df) <- c("event_id", "gene_id","gene", "chr", "event_type", "FDR", "IncLevel_CNT", "IncLevel_PMA", "IncLevelDiff")
MXE_df <- dplyr::select(MXE, ID...1, GeneID, geneSymbol, chr, ID...12, FDR, IncLevel1, IncLevel2, IncLevelDifference)
colnames(MXE_df) <- c("event_id","gene_id", "gene", "chr", "event_type", "FDR", "IncLevel_CNT", "IncLevel_PMA", "IncLevelDiff")
RI_df <- dplyr::select(RI, ID...1, GeneID, geneSymbol, chr, ID...12, FDR, IncLevel1, IncLevel2, IncLevelDifference)
colnames(RI_df) <- c("event_id", "gene_id", "gene", "chr", "event_type", "FDR", "IncLevel_CNT", "IncLevel_PMA", "IncLevelDiff")
A3SS_df <- dplyr::select(A3SS, ID...1, GeneID, geneSymbol, chr, ID...12, FDR, IncLevel1, IncLevel2, IncLevelDifference)
colnames(A3SS_df) <- c("event_id", "gene_id", "gene", "chr", "event_type", "FDR", "IncLevel_CNT", "IncLevel_PMA", "IncLevelDiff")
A5SS_df <- dplyr::select(A5SS, ID...1, GeneID, geneSymbol, chr, ID...12, FDR, IncLevel1, IncLevel2, IncLevelDifference)
colnames(A5SS_df) <- c("event_id", "gene_id", "gene", "chr", "event_type", "FDR", "IncLevel_CNT", "IncLevel_PMA", "IncLevelDiff")

rmats <- rbind(SE_df, MXE_df, RI_df, A3SS_df, A5SS_df) # merge tables
saveRDS(rmats, "DASE_results.rds") # save rds for future use

rmats_df <- rmats %>%
  mutate(IJC_SAMPLE_1 = strsplit(as.character(IJC_SAMPLE_1), ","),
         SJC_SAMPLE_1 = strsplit(as.character(SJC_SAMPLE_1), ","),
         IJC_SAMPLE_2 = strsplit(as.character(IJC_SAMPLE_2), ","),
         SJC_SAMPLE_2 = strsplit(as.character(SJC_SAMPLE_2), ",")) %>%
  unnest_wider(IJC_SAMPLE_1, names_sep = "_rep") %>%
  unnest_wider(SJC_SAMPLE_1, names_sep = "_rep") %>%
  unnest_wider(IJC_SAMPLE_2, names_sep = "_rep") %>%
  unnest_wider(SJC_SAMPLE_2, names_sep = "_rep") # separate counts for each replicate

rmats_df <- rmats_df %>% mutate(across(matches("rep"), as.numeric)) # convert to numerical values

rmats_df <- rmats_df %>%
  mutate(across(matches("IJC_SAMPLE_1_rep|SJC_SAMPLE_1_rep"), ~replace_na(.,0))) %>%
  mutate(across(matches("IJC_SAMPLE_2_rep|SJC_SAMPLE_2_rep"), ~replace_na(.,0))) # get total counts per replicate

rmats_df <- rmats_df %>%
  rowwise() %>% mutate(TotalCounts_mono = mean(c_across(matches("IJC_SAMPLE_1_rep|SJC_SAMPLE_1_rep"))),
                       TotalCounts_macro = mean(c_across(matches("IJC_SAMPLE_2_rep|SJC_SAMPLE_2_rep")))) %>% ungroup()

rmats_df$MeanTotal <- rowMeans(rmats_df[,c("TotalCounts_mono","TotalCounts_macro")])

rmats_filter <- rmats_df %>% filter(MeanTotal >= 10) # filter out reads with low average total counts (rats default)

rmats_counts_long <- rmats_filter %>%
  select(TotalCounts_mono, TotalCounts_macro) %>%
  mutate(EventID = row_number()) %>%
  melt(id.vars="EventID", variable.name="sample", value.name="counts")

rmats_counts_long$condition <- factor(ifelse(rmats_counts_long$sample == "TotalCounts_mono", "Monocyte", "Macrophage"),
                                      levels=c("Monocyte","Macrophage"))

rmats_counts_long$logCounts <- log2(rmats_counts_long$counts + 1)

das_counts_violin <- ggplot(rmats_counts_long, aes(x = condition, y = logCounts, fill = condition)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
  stat_summary(fun = "median", geom = "point", size = 2, color = "black") +
  labs(x = "Cell Type", y = "log2(Junction Counts)") +
  scale_fill_manual(values = c("Monocyte" = "#46B1E1", "Macrophage" = "#78206E")) +
  plot_theme() # making violin plots to show expression bias 

ggsave("das_junccounts_violin.png", das_counts_violin, units = "mm", height = 80, width = 100 , dpi = 500)

#### Filtering of aligned reads ####
bam_outcomes <- read_tsv("CNTvsPMA/2024-07-03-14_15_47_619298_read_outcomes_by_bam.txt", col_names = FALSE)

# divide file by replicate
bam1 <- bam_outcomes %>% dplyr::slice(1:11)
bam1 <- separate_wider_delim(bam1, cols = X1, delim = ": ", names_sep = "", too_few = "align_start")
CNT1 <- pivot_wider(bam1, names_from = X11, values_from = X12, names_prefix = "")
bam2 <- bam_outcomes %>% dplyr::slice(12:22)
bam2 <- separate_wider_delim(bam2, cols = X1, delim = ": ", names_sep = "", too_few = "align_start")
CNT2 <- pivot_wider(bam2, names_from = X11, values_from = X12, names_prefix = "")
bam3 <- bam_outcomes %>% dplyr::slice(23:33)
bam3 <- separate_wider_delim(bam3, cols = X1, delim = ": ", names_sep = "", too_few = "align_start")
CNT3 <- pivot_wider(bam3, names_from = X11, values_from = X12, names_prefix = "")
bam4 <- bam_outcomes %>% dplyr::slice(34:44)
bam4 <- separate_wider_delim(bam4, cols = X1, delim = ": ", names_sep = "", too_few = "align_start")
PMA1 <- pivot_wider(bam4, names_from = X11, values_from = X12, names_prefix = "")
bam5 <- bam_outcomes %>% dplyr::slice(45:55)
bam5 <- separate_wider_delim(bam5, cols = X1, delim = ": ", names_sep = "", too_few = "align_start")
PMA2 <- pivot_wider(bam5, names_from = X11, values_from = X12, names_prefix = "")
bam6 <- bam_outcomes %>% dplyr::slice(56:66)
bam6 <- separate_wider_delim(bam6, cols = X1, delim = ": ", names_sep = "", too_few = "align_start")
PMA3 <- pivot_wider(bam6, names_from = X11, values_from = X12, names_prefix = "")

# extract BAM counts and assign sample name
CNT1_counts <- dplyr::select(CNT1, USED, TOTAL_FOR_BAM) 
CNT1_counts$sample <- "CNT1"
CNT2_counts <- dplyr::select(CNT2, USED, TOTAL_FOR_BAM)
CNT2_counts$sample <- "CNT2"
CNT3_counts <- dplyr::select(CNT3, USED, TOTAL_FOR_BAM)
CNT3_counts$sample <- "CNT3"

PMA1_counts <- dplyr::select(PMA1, USED, TOTAL_FOR_BAM)
PMA1_counts$sample <- "PMA1"
PMA2_counts <- dplyr::select(PMA2, USED, TOTAL_FOR_BAM)
PMA2_counts$sample <- "PMA2"
PMA3_counts <- dplyr::select(PMA3, USED, TOTAL_FOR_BAM)
PMA3_counts$sample <- "PMA3"

BAM_counts <- rbind(CNT1_counts, CNT2_counts, CNT3_counts, PMA1_counts, PMA2_counts, PMA3_counts) # combine samples into single df for plotting 
BAM_counts_long <- pivot_longer(BAM_counts, cols = c(USED, TOTAL_FOR_BAM), names_to = "Count")

bam_outcome <- ggplot(BAM_counts_long, aes(x = sample, y = value, fill = Count)) + geom_bar(stat = "identity",position = "dodge") +
  scale_fill_manual(values = c("cadetblue", "blueviolet"), labels = c("Total Reads", "Used Reads")) + 
  labs(y = "Read counts", x = "sample", fill = "BAM outcome") + 
  plot_theme()

ggsave("BAM_outcomes.png", bam_outcome,  units = "mm", height = 200 , width = 250, dpi = 500)


### filtering criteria ###
