# R script for DTU analysis using DRIMSeq and quality control metrics (transcript abundance distribution and transcriptome complexity)

# load libraries 
library(tidyverse)
library(DRIMSeq)
library(tximport)
library(GenomicFeatures)
library(txdbmaker)
library(stageR)

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

#### DTU analysis ####
sampleTable <- data.frame(sample_id = c("CNT1","CNT2","CNT3","PMA1", "PMA2","PMA3"),
                          condition = as.factor(c("Control", "Control","Control","PMA", "PMA","PMA")),
                          libType = "paired-end") # create Sample table with sample information

files <- c("RSEM/CNT_replicate1/CNT_replicate_1.isoforms.results",
           "RSEM/CNT_replicate2/CNT_replicate_2.isoforms.results",
           "RSEM/CNT_replicate3/CNT_replicate_3.isoforms.results",
           "RSEM/PMA_replicate1/PMA_rep1.isoforms.results",
           "RSEM/PMA_replicate2/PMA_rep2.isoforms.results",
           "RSEM/PMA_replicate3/PMA_rep3.isoforms.results")
names(files) <- sampleTable$sample_id # create sample table with estimated transcript abundances

# extract gene and transcript names from gtf file used for quantification
txdb <- makeTxDbFromGFF("gencode.v44.primary_assembly.annotation.gtf/gencode.v44.primary_assembly.annotation.gtf")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
counts <- tximport(files, type = "rsem", txIn = TRUE ,txOut = TRUE,
                   tx2gene = tx2gene, countsFromAbundance = "dtuScaledTPM") # collapse transcripts to gene level and scale for DTU analysis

count <- counts$counts # extract count matrix from tximport matrix
countTable <- count[rowSums(count) > 0, ] #filter non-expressed transcripts
tx2gene_sub <- tx2gene[match(rownames(countTable),tx2gene$TXNAME), ] # create count matrix to create DRIMSeq data set
countTable <- data.frame(gene_id = tx2gene_sub$GENEID, feature_id = tx2gene_sub$TXNAME,countTable)
dmDS <- dmDSdata(counts = countTable, samples = sampleTable) # create dmDSdata object for DRIMSeq analysis 
#filter lowly expressed transcripts
n = nrow(sampleTable)
n.small = min(table(sampleTable$condition))
dmDS <- dmFilter(dmDS,
                 min_samps_feature_expr = n.small, min_feature_expr = 10,
                 min_samps_feature_prop = n.small, min_feature_prop = 0.1,
                 min_samps_gene_expr = n, min_gene_expr = 10) 

design <- model.matrix(~condition, data = DRIMSeq::samples(dmDS)) # create design matrix using sample data and design formula 

dmDS <- dmPrecision(dmDS, design = design) # calculate precision estimates 
dmDS <- dmFit(dmDS, design = design) # fit regression coefficients, gene-level estimates and transcript-level bb model results
dmDS <- dmTest(dmDS, coef = "conditionPMA") # perform null hypothesis testing between conditions
# extract genes that show evidence of DTU and the differentially used transcripts
results_gene <- DRIMSeq::results(dmDS)
results_transcripts <- DRIMSeq::results(dmDS, level = "feature")

#### Stage wise analysis ####
pval_screen <- results_gene$pvalue # extract p-values of genes exhibiting DTU
names(pval_screen) <- results_gene$gene_id
pval_screen <- ifelse(is.na(pval_screen),1,pval_screen) # remove NA values 
pval_confirmation <- matrix(results_transcripts$pvalue, ncol = 1) # extract p-values of differentially used transcripts
rownames(pval_confirmation) <- results_transcripts$feature_id

pval_confirmation <- ifelse(is.na(pval_confirmation),1,pval_confirmation) # remove NA values
txs2gene <- results_transcripts[ , c("feature_id","gene_id")]
stage_Object <- stageRTx(pval_screen, pval_confirmation, pScreenAdjusted = FALSE, tx2gene = txs2gene)
stage_Object <- stageWiseAdjustment(stage_Object, method = 'dtu', alpha = 0.05) # perform test at OFDR of 0.05
padj <- getAdjustedPValues(stage_Object, order = FALSE, onlySignificantGenes = FALSE)

saveRDS(padj, 'DTU_stageR_results.rds') # save RDS for downstream analysis

#### estimated transcript abundance distribution ####
# extract transcript abundances and convert into log values
tpm <- counts$abundance 
logtpm <- log2(tpm + 1) # shift + 1

TPM_long <- melt(logtpm, varnames=c("transcript","sample"), value.name="logTPM")
TPM_long$condition <- factor(ifelse(grepl("CNT", TPM_long$sample), "Monocyte", "Macrophage"), levels = c("Monocyte", "Macrophage"))

dtu_tpm_violin <- ggplot(TPM_long, aes(x = condition, y = logTPM, fill = condition)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
  stat_summary(fun = "median", geom = "point", size = 2, color = "black") +
  labs(x = "Cell Type", y = "log2(TPM)") +
  scale_fill_manual(values = c("Monocyte" = "#46B1E1", "Macrophage" = "#78206E")) +
  plot_theme() ## making violin plots to show expression bias 

ggsave("trasncript_distribution.png", dtu_tpm_violin, units = "mm", height = 80, width = 100 , dpi = 500)

#### transcriptome complexity ####
prop <- proportions(dmDS) # obtain transcript 

txs2gene <- data.frame(feature_id = prop$feature_id,
                       gene_id = prop$gene_id) # match transcripts 2 genes

var_per_tx <- prop %>%
  select(feature_id, gene_id, matches("CNT|PMA")) %>%
  rowwise() %>%
  mutate(var_tx = var(c_across(matches("CNT|PMA")))) %>%
  ungroup() # Calculate variance of proportions per transcript across samples

var_df <- var_per_tx %>%
  group_by(gene_id) %>%
  summarise(n_tx = n_distinct(feature_id),
            var_props = mean(var_tx)) # summarise to gene level

transcript_complex <- ggplot(var_df, aes(x = n_tx, y = var_props)) +
  geom_point(alpha=0.5) +
  geom_smooth(method="loess") +
  labs(x="Number of transcripts per gene",
       y="Mean variance of transcript proportions") +
  plot_theme()
ggsave("Transcriptome_complexity.png", transcript_complex, units = "mm", height = 100, width = 130 , dpi = 500)

