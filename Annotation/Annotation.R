setwd("C:/Users/huber/Box/Dave and Ehren Collaborations and stuff/Tsp_Thermal_Paper/Annotation")

#clear environment
rm(list=ls())+
  dev.off()

# load annotation file
info <- read.csv("Tsp_transcript_info_updated.csv", header=TRUE)
# load DE results 
DE <- read.delim("All_DE_Transcripts_padj0.05.tsv", header = TRUE, sep = "\t")

# create new column DE$gene that pulls the gene name from info by transcript ID
DE$gene <- info$gene[match(DE$ID, info$ID)]

# clean up gene names
DE$gene <- gsub("_", " ", DE$gene)

# what % of transcripts have annotation?
library(dplyr)

#### write annotated DE results to file ####

write.csv(DE, file = "DE_Results_All_annotated.csv", row.names = FALSE)
