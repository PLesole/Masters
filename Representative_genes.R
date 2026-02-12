# R script to generate figures for representative genes in DEXSeq and DRIMSeq

# load libraries 
library(DEXSeq)
library(DRIMSeq)
library(tidyverse)

# load output files 
DEU <- readRDS('DEU_results.rds')
DTU <- readRDS('DTU_dmDS.rds')

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

#### plot DEU ####
DEXSeq::plotDEXSeq(DEU, "ENSG00000133048.13", names = TRUE, legend = TRUE, splicing = TRUE, color = c("#46B1E1","#78206E")) # CHI3L1
DEXSeq::plotDEXSeq(DEU, "ENSG00000130203.10", names = TRUE, legend = TRUE, splicing = TRUE, color = c("#46B1E1","#78206E")) # APOE
DEXSeq::plotDEXSeq(DEU, "ENSG00000269404.7", names = TRUE, legend = TRUE, splicing = TRUE, color = c("#46B1E1","#78206E")) # SPIB


#### plot DTU ####
proportions <- DRIMSeq::proportions(DTU) 
###### APOE ####
proportions_filter <- filter(proportions, gene_id == "ENSG00000130203.10")
proportions_filter$transcript_id <- c("APOE_201", "APOE_203") 
proportions_filter <- mutate(proportions_filter, CNT = (CNT1 + CNT2 + CNT3)/3,
                             PMA = (PMA1 + PMA2 + PMA3)/3) # average proportions per cell type
proportions_filter_long <- pivot_longer(proportions_filter, cols = c("CNT","PMA") , names_to = "Condition" , values_to = "ave_proportion")

apoe_proportion <- ggplot(proportions_filter_long, aes(x = transcript_id, y = ave_proportion , fill = Condition)) + geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#46B1E1","#78206E"), labels = c("Monocytes","Macrophages")) + 
  scale_y_continuous(expand = c(0,0)) +  
  labs(x = "Transcript", y = "Average Relative Proportion", fill = "Cell Type") + 
  plot_theme()
ggsave("DTU_APOE_proportions.png", apoe_proportion, units = "mm", height = 90, width = 120 ,dpi = 500)

###### CHI3L1- unable to compute ####  
#proportions_filter <- filter(proportions, gene_id == "ENSG00000133048.13")
#proportions_filter$transcript_id <- c("", "") 
#proportions_filter <- mutate(proportions_filter, CNT = (CNT1 + CNT2 + CNT3)/3,
#                             PMA = (PMA1 + PMA2 + PMA3)/3) # average proportions per cell type
#proportions_filter_long <- pivot_longer(proportions_filter, cols = c("CNT","PMA") , names_to = "Condition" , values_to = "ave_proportion")

#chi3l1_proportion <- ggplot(proportions_filter_long, aes(x = transcript_id, y = ave_proportion , fill = Condition)) + geom_bar(stat = "identity", position = "dodge") +
#  scale_fill_manual(values = c("#46B1E1","#78206E"), labels = c("Monocytes","Macrophages")) + 
#  scale_y_continuous(expand = c(0,0)) +  
#  labs(x = "Transcript", y = "Average Relative Proportion", fill = "Cell Type") + 
#  plot_theme()
#ggsave("DTU_CHI3L1_proportions.png", chi3l1_proportion, units = "mm", height = 80, width = 100 ,dpi = 500)

##### SPIB ####
proportions_filter <- filter(proportions, gene_id == "ENSG00000269404.7")
proportions_filter$transcript_id <- c("SPIB_201", "SPIB_206","SPIB_207") 
proportions_filter <- mutate(proportions_filter, CNT = (CNT1 + CNT2 + CNT3)/3,
                             PMA = (PMA1 + PMA2 + PMA3)/3) # average proportions per cell type
proportions_filter_long <- pivot_longer(proportions_filter, cols = c("CNT","PMA") , names_to = "Condition" , values_to = "ave_proportion")

spib_proportion <- ggplot(proportions_filter_long, aes(x = transcript_id, y = ave_proportion , fill = Condition)) + geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#46B1E1","#78206E"), labels = c("Monocytes","Macrophages")) + 
  scale_y_continuous(expand = c(0,0)) +  
  labs(x = "Transcript", y = "Average Relative Proportion", fill = "Cell Type") + 
  plot_theme()
ggsave("DTU_SPIB_proportions.png", spib_proportion, units = "mm", height = 90, width = 120 ,dpi = 500)
