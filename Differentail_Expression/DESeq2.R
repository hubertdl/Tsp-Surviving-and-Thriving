{
#Set working directory
setwd("")

# Load libraries
library(DESeq2)
library(heatmap3)
library(EnhancedVolcano)
library(RColorBrewer)
library(ggrepel)
library(apeglm)
library(ashr)

  # Clear environment
  rm(list = ls())
  graphics.off()
  cat("\014")
  
  # Load the full counts.tab file
  df <- read.table("Expression_tables/Tsp_thermal_counts.tsv", header = TRUE, sep = "\t", row.names = 1)
  head(df)
  
  # Load the full key.tab file
  dkey <- read.table("Expression_tables/Tsp_thermal_KEY.tsv", header=T, sep = "\t", stringsAsFactors = FALSE)
  dkey$CONDITION <- as.factor(dkey$CONDITION)
  dkey$TISSUE <- as.factor(dkey$TISSUE)
  dkey
  
  # Create output directories
  dir.create("Heatmaps/Wald_Heatmaps", showWarnings = FALSE, recursive = TRUE)
  dir.create("Results/Wald_Results", showWarnings = FALSE, recursive = TRUE)
  dir.create("Volcano_Plots", showWarnings = FALSE)
  dir.create("Results/Wald_Results/ShrunkenL2FC", showWarnings = FALSE, recursive = TRUE)
  
  # Set thresholds and colors
  coveragethrd <- 3 # min average number of reads per sample 
  WALD_FDR = 0.05
  WALD_lfc_cutoff = .585
  hmcol<-colorRampPalette(c("Blue","Black","Red"))(15) # color scale for heatmaps
}

############################################
# Subset data by TISSUE 
############################################
{
#  columns <- ncol(df)
#  df <- df[rowSums(df)>=(columns*coveragethrd),]
#  print(paste("Genes after filtering:", nrow(df)))
  
  # BRAIN - COLD vs CONTROL
  dkey_BRAIN_COLD <- dkey[dkey$TISSUE == "BRAIN" & dkey$CONDITION %in% c("COLD", "CONTROL"), ]
  df_BRAIN_COLD <- df[, dkey_BRAIN_COLD$ID]
  columns_BRAIN_COLD <- ncol(df_BRAIN_COLD)
  df_BRAIN_COLD <- df_BRAIN_COLD[rowSums(df_BRAIN_COLD)>=(columns_BRAIN_COLD*coveragethrd),]
  print(paste("BRAIN COLD vs CONTROL genes after filtering:", nrow(df_BRAIN_COLD)))
  print(table(dkey_BRAIN_COLD$CONDITION))
  
  # BRAIN - HEAT vs CONTROL
  dkey_BRAIN_HEAT <- dkey[dkey$TISSUE == "BRAIN" & dkey$CONDITION %in% c("HEAT", "CONTROL"), ]
  df_BRAIN_HEAT <- df[, dkey_BRAIN_HEAT$ID]
  columns_BRAIN_HEAT <- ncol(df_BRAIN_HEAT)
  df_BRAIN_HEAT <- df_BRAIN_HEAT[rowSums(df_BRAIN_HEAT)>=(columns_BRAIN_HEAT*coveragethrd),]
  print(paste("BRAIN HEAT vs CONTROL genes after filtering:", nrow(df_BRAIN_HEAT)))
  print(table(dkey_BRAIN_HEAT$CONDITION))
  
  # HEART - COLD vs CONTROL
  dkey_HEART_COLD <- dkey[dkey$TISSUE == "HEART" & dkey$CONDITION %in% c("COLD", "CONTROL"), ]
  df_HEART_COLD <- df[, dkey_HEART_COLD$ID]
  columns_HEART_COLD <- ncol(df_HEART_COLD)
  df_HEART_COLD <- df_HEART_COLD[rowSums(df_HEART_COLD)>=(columns_HEART_COLD*coveragethrd),]
  print(paste("HEART COLD vs CONTROL genes after filtering:", nrow(df_HEART_COLD)))
  print(table(dkey_HEART_COLD$CONDITION))
  
  # HEART - HEAT vs CONTROL
  dkey_HEART_HEAT <- dkey[dkey$TISSUE == "HEART" & dkey$CONDITION %in% c("HEAT", "CONTROL"), ]
  df_HEART_HEAT <- df[, dkey_HEART_HEAT$ID]
  columns_HEART_HEAT <- ncol(df_HEART_HEAT)
  df_HEART_HEAT <- df_HEART_HEAT[rowSums(df_HEART_HEAT)>=(columns_HEART_HEAT*coveragethrd),]
  print(paste("HEART HEAT vs CONTROL genes after filtering:", nrow(df_HEART_HEAT)))
  print(table(dkey_HEART_HEAT$CONDITION))
  
  # LIVER - COLD vs CONTROL
  dkey_LIVER_COLD <- dkey[dkey$TISSUE == "LIVER" & dkey$CONDITION %in% c("COLD", "CONTROL"), ]
  df_LIVER_COLD <- df[, dkey_LIVER_COLD$ID]
  columns_LIVER_COLD <- ncol(df_LIVER_COLD)
  df_LIVER_COLD <- df_LIVER_COLD[rowSums(df_LIVER_COLD)>=(columns_LIVER_COLD*coveragethrd),]
  print(paste("LIVER COLD vs CONTROL genes after filtering:", nrow(df_LIVER_COLD)))
  print(table(dkey_LIVER_COLD$CONDITION))
  
  # LIVER - HEAT vs CONTROL
  dkey_LIVER_HEAT <- dkey[dkey$TISSUE == "LIVER" & dkey$CONDITION %in% c("HEAT", "CONTROL"), ]
  df_LIVER_HEAT <- df[, dkey_LIVER_HEAT$ID]
  columns_LIVER_HEAT <- ncol(df_LIVER_HEAT)
  df_LIVER_HEAT <- df_LIVER_HEAT[rowSums(df_LIVER_HEAT)>=(columns_LIVER_HEAT*coveragethrd),]
  print(paste("LIVER HEAT vs CONTROL genes after filtering:", nrow(df_LIVER_HEAT)))
  print(table(dkey_LIVER_HEAT$CONDITION))
  
  # TESTIS - COLD vs CONTROL
  dkey_TESTIS_COLD <- dkey[dkey$TISSUE == "TESTIS" & dkey$CONDITION %in% c("COLD", "CONTROL"), ]
  df_TESTIS_COLD <- df[, dkey_TESTIS_COLD$ID]
  columns_TESTIS_COLD <- ncol(df_TESTIS_COLD)
  df_TESTIS_COLD <- df_TESTIS_COLD[rowSums(df_TESTIS_COLD)>=(columns_TESTIS_COLD*coveragethrd),]
  print(paste("TESTIS COLD vs CONTROL genes after filtering:", nrow(df_TESTIS_COLD)))
  print(table(dkey_TESTIS_COLD$CONDITION))
  
  # TESTIS - HEAT vs CONTROL
  dkey_TESTIS_HEAT <- dkey[dkey$TISSUE == "TESTIS" & dkey$CONDITION %in% c("HEAT", "CONTROL"), ]
  df_TESTIS_HEAT <- df[, dkey_TESTIS_HEAT$ID]
  columns_TESTIS_HEAT <- ncol(df_TESTIS_HEAT)
  df_TESTIS_HEAT <- df_TESTIS_HEAT[rowSums(df_TESTIS_HEAT)>=(columns_TESTIS_HEAT*coveragethrd),]
  print(paste("TESTIS HEAT vs CONTROL genes after filtering:", nrow(df_TESTIS_HEAT)))
  print(table(dkey_TESTIS_HEAT$CONDITION))
}
  
################################################################################
# WALD TESTS FOR ALL PAIRWISE COMPARISONS
################################################################################
{
  #######################
  # WALD test for COLD vs CONTROL in BRAIN
  #######################
  dm_BRAIN_COLD <- as.matrix(df_BRAIN_COLD)
  dds_BRAIN_COLD <- DESeqDataSetFromMatrix(dm_BRAIN_COLD, dkey_BRAIN_COLD, design = ~ CONDITION)
  dds_BRAIN_COLD$CONDITION <- relevel(dds_BRAIN_COLD$CONDITION, ref="CONTROL")
  dds_BRAIN_COLD <- DESeq(dds_BRAIN_COLD, test = "Wald", parallel = TRUE)
  res_BRAIN_COLD <- results(dds_BRAIN_COLD, test = "Wald", alpha = WALD_FDR)
  res_BRAIN_COLD <- res_BRAIN_COLD[order(res_BRAIN_COLD$padj), ]
  nrow(res_BRAIN_COLD[which(res_BRAIN_COLD$padj <= WALD_FDR),])
  summary(res_BRAIN_COLD)
  
  #######################
  # WALD test for HEAT vs CONTROL in BRAIN
  #######################
  dm_BRAIN_HEAT <- as.matrix(df_BRAIN_HEAT)
  dds_BRAIN_HEAT <- DESeqDataSetFromMatrix(dm_BRAIN_HEAT, dkey_BRAIN_HEAT, design = ~ CONDITION)
  dds_BRAIN_HEAT$CONDITION <- relevel(dds_BRAIN_HEAT$CONDITION, ref="CONTROL")
  dds_BRAIN_HEAT <- DESeq(dds_BRAIN_HEAT, test = "Wald", parallel = TRUE)
  res_BRAIN_HEAT <- results(dds_BRAIN_HEAT, test = "Wald", alpha = WALD_FDR)
  res_BRAIN_HEAT <- res_BRAIN_HEAT[order(res_BRAIN_HEAT$padj), ]
  nrow(res_BRAIN_HEAT[which(res_BRAIN_HEAT$padj <= WALD_FDR),])
  summary(res_BRAIN_HEAT)
  
  #######################
  # WALD test for COLD vs CONTROL in HEART
  #######################
  dm_HEART_COLD <- as.matrix(df_HEART_COLD)
  dds_HEART_COLD <- DESeqDataSetFromMatrix(dm_HEART_COLD, dkey_HEART_COLD, design = ~ CONDITION)
  dds_HEART_COLD$CONDITION <- relevel(dds_HEART_COLD$CONDITION, ref="CONTROL")
  dds_HEART_COLD <- DESeq(dds_HEART_COLD, test = "Wald", parallel = TRUE)
  res_HEART_COLD <- results(dds_HEART_COLD, test = "Wald", alpha = WALD_FDR)
  res_HEART_COLD <- res_HEART_COLD[order(res_HEART_COLD$padj), ]
  nrow(res_HEART_COLD[which(res_HEART_COLD$padj <= WALD_FDR),])
  summary(res_HEART_COLD)
  
  #######################
  # WALD test for HEAT vs CONTROL in HEART
  #######################
  dm_HEART_HEAT <- as.matrix(df_HEART_HEAT)
  dds_HEART_HEAT <- DESeqDataSetFromMatrix(dm_HEART_HEAT, dkey_HEART_HEAT, design = ~ CONDITION)
  dds_HEART_HEAT$CONDITION <- relevel(dds_HEART_HEAT$CONDITION, ref="CONTROL")
  dds_HEART_HEAT <- DESeq(dds_HEART_HEAT, test = "Wald", parallel = TRUE)
  res_HEART_HEAT <- results(dds_HEART_HEAT, test = "Wald", alpha = WALD_FDR)
  res_HEART_HEAT <- res_HEART_HEAT[order(res_HEART_HEAT$padj), ]
  nrow(res_HEART_HEAT[which(res_HEART_HEAT$padj <= WALD_FDR),])
  summary(res_HEART_HEAT)
  
  #######################
  # WALD test for COLD vs CONTROL in LIVER
  #######################
  dm_LIVER_COLD <- as.matrix(df_LIVER_COLD)
  dds_LIVER_COLD <- DESeqDataSetFromMatrix(dm_LIVER_COLD, dkey_LIVER_COLD, design = ~ CONDITION)
  dds_LIVER_COLD$CONDITION <- relevel(dds_LIVER_COLD$CONDITION, ref="CONTROL")
  dds_LIVER_COLD <- DESeq(dds_LIVER_COLD, test = "Wald", parallel = TRUE)
  res_LIVER_COLD <- results(dds_LIVER_COLD, test = "Wald", alpha = WALD_FDR)
  res_LIVER_COLD <- res_LIVER_COLD[order(res_LIVER_COLD$padj), ]
  nrow(res_LIVER_COLD[which(res_LIVER_COLD$padj <= WALD_FDR),])
  summary(res_LIVER_COLD)
  
  #######################
  # WALD test for HEAT vs CONTROL in LIVER
  #######################
  dm_LIVER_HEAT <- as.matrix(df_LIVER_HEAT)
  dds_LIVER_HEAT <- DESeqDataSetFromMatrix(dm_LIVER_HEAT, dkey_LIVER_HEAT, design = ~ CONDITION)
  dds_LIVER_HEAT$CONDITION <- relevel(dds_LIVER_HEAT$CONDITION, ref="CONTROL")
  dds_LIVER_HEAT <- DESeq(dds_LIVER_HEAT, test = "Wald", parallel = TRUE)
  res_LIVER_HEAT <- results(dds_LIVER_HEAT, test = "Wald", alpha = WALD_FDR)
  res_LIVER_HEAT <- res_LIVER_HEAT[order(res_LIVER_HEAT$padj), ]
  nrow(res_LIVER_HEAT[which(res_LIVER_HEAT$padj <= WALD_FDR),])
  summary(res_LIVER_HEAT)
  
  #######################
  # WALD test for COLD vs CONTROL in TESTIS
  #######################
  dm_TESTIS_COLD <- as.matrix(df_TESTIS_COLD)
  dds_TESTIS_COLD <- DESeqDataSetFromMatrix(dm_TESTIS_COLD, dkey_TESTIS_COLD, design = ~ CONDITION)
  dds_TESTIS_COLD$CONDITION <- relevel(dds_TESTIS_COLD$CONDITION, ref="CONTROL")
  dds_TESTIS_COLD <- DESeq(dds_TESTIS_COLD, test = "Wald", parallel = TRUE)
  res_TESTIS_COLD <- results(dds_TESTIS_COLD, test = "Wald", alpha = WALD_FDR)
  res_TESTIS_COLD <- res_TESTIS_COLD[order(res_TESTIS_COLD$padj), ]
  nrow(res_TESTIS_COLD[which(res_TESTIS_COLD$padj <= WALD_FDR),])
  summary(res_TESTIS_COLD)
  
  #######################
  # WALD test for HEAT vs CONTROL in TESTIS
  #######################
  dm_TESTIS_HEAT <- as.matrix(df_TESTIS_HEAT)
  dds_TESTIS_HEAT <- DESeqDataSetFromMatrix(dm_TESTIS_HEAT, dkey_TESTIS_HEAT, design = ~ CONDITION)
  dds_TESTIS_HEAT$CONDITION <- relevel(dds_TESTIS_HEAT$CONDITION, ref="CONTROL")
  dds_TESTIS_HEAT <- DESeq(dds_TESTIS_HEAT, test = "Wald", parallel = TRUE)
  res_TESTIS_HEAT <- results(dds_TESTIS_HEAT, test = "Wald", alpha = WALD_FDR)
  res_TESTIS_HEAT <- res_TESTIS_HEAT[order(res_TESTIS_HEAT$padj), ]
  nrow(res_TESTIS_HEAT[which(res_TESTIS_HEAT$padj <= WALD_FDR),])
  summary(res_TESTIS_HEAT)
}

#######################
# Export Results
#######################
{  
  # BRAIN_COLD
  if(exists("res_BRAIN_COLD") && !is.null(res_BRAIN_COLD)) {
    print(paste("DEGs at FDR ≤", WALD_FDR, ":", nrow(res_BRAIN_COLD[which(res_BRAIN_COLD$padj <= WALD_FDR),])))
    print(paste("DEGs at FDR ≤", WALD_FDR, "& |LFC| ≥", WALD_lfc_cutoff, ":", nrow(res_BRAIN_COLD[which(res_BRAIN_COLD$padj <= WALD_FDR & abs(res_BRAIN_COLD$log2FoldChange) >= WALD_lfc_cutoff),])))
    summary(res_BRAIN_COLD)
    write.table(res_BRAIN_COLD, "Results/Wald_Results/BRAIN_COLD_vs_CONTROL_WALD_DE.tsv", sep="\t", col.names=NA, quote = FALSE)
  } else {
    print("res_BRAIN_COLD not found or is NULL")
  }
  
  # BRAIN_HEAT
  if(exists("res_BRAIN_HEAT") && !is.null(res_BRAIN_HEAT)) {
    print(paste("DEGs at FDR ≤", WALD_FDR, ":", nrow(res_BRAIN_HEAT[which(res_BRAIN_HEAT$padj <= WALD_FDR),])))
    print(paste("DEGs at FDR ≤", WALD_FDR, "& |LFC| ≥", WALD_lfc_cutoff, ":", nrow(res_BRAIN_HEAT[which(res_BRAIN_HEAT$padj <= WALD_FDR & abs(res_BRAIN_HEAT$log2FoldChange) >= WALD_lfc_cutoff),])))
    summary(res_BRAIN_HEAT)
    write.table(res_BRAIN_HEAT, "Results/Wald_Results/BRAIN_HEAT_vs_CONTROL_WALD_DE.tsv", sep="\t", col.names=NA, quote = FALSE)
  } else {
    print("res_BRAIN_HEAT not found or is NULL")
  }
  
  # HEART_COLD
  if(exists("res_HEART_COLD") && !is.null(res_HEART_COLD)) {
    print(paste("DEGs at FDR ≤", WALD_FDR, ":", nrow(res_HEART_COLD[which(res_HEART_COLD$padj <= WALD_FDR),])))
    print(paste("DEGs at FDR ≤", WALD_FDR, "& |LFC| ≥", WALD_lfc_cutoff, ":", nrow(res_HEART_COLD[which(res_HEART_COLD$padj <= WALD_FDR & abs(res_HEART_COLD$log2FoldChange) >= WALD_lfc_cutoff),])))
    summary(res_HEART_COLD)
    write.table(res_HEART_COLD, "Results/Wald_Results/HEART_COLD_vs_CONTROL_WALD_DE.tsv", sep="\t", col.names=NA, quote = FALSE)
  } else {
    print("res_HEART_COLD not found or is NULL")
  }
  
  # HEART_HEAT
  if(exists("res_HEART_HEAT") && !is.null(res_HEART_HEAT)) {
    print(paste("DEGs at FDR ≤", WALD_FDR, ":", nrow(res_HEART_HEAT[which(res_HEART_HEAT$padj <= WALD_FDR),])))
    print(paste("DEGs at FDR ≤", WALD_FDR, "& |LFC| ≥", WALD_lfc_cutoff, ":", nrow(res_HEART_HEAT[which(res_HEART_HEAT$padj <= WALD_FDR & abs(res_HEART_HEAT$log2FoldChange) >= WALD_lfc_cutoff),])))
    summary(res_HEART_HEAT)
    write.table(res_HEART_HEAT, "Results/Wald_Results/HEART_HEAT_vs_CONTROL_WALD_DE.tsv", sep="\t", col.names=NA, quote = FALSE)
  } else {
    print("res_HEART_HEAT not found or is NULL")
  }
  
  # LIVER_COLD
  if(exists("res_LIVER_COLD") && !is.null(res_LIVER_COLD)) {
    print(paste("DEGs at FDR ≤", WALD_FDR, ":", nrow(res_LIVER_COLD[which(res_LIVER_COLD$padj <= WALD_FDR),])))
    print(paste("DEGs at FDR ≤", WALD_FDR, "& |LFC| ≥", WALD_lfc_cutoff, ":", nrow(res_LIVER_COLD[which(res_LIVER_COLD$padj <= WALD_FDR & abs(res_LIVER_COLD$log2FoldChange) >= WALD_lfc_cutoff),])))
    summary(res_LIVER_COLD)
    write.table(res_LIVER_COLD, "Results/Wald_Results/LIVER_COLD_vs_CONTROL_WALD_DE.tsv", sep="\t", col.names=NA, quote = FALSE)
  } else {
    print("res_LIVER_COLD not found or is NULL")
  }
  
  # LIVER_HEAT
  if(exists("res_LIVER_HEAT") && !is.null(res_LIVER_HEAT)) {
    print(paste("DEGs at FDR ≤", WALD_FDR, ":", nrow(res_LIVER_HEAT[which(res_LIVER_HEAT$padj <= WALD_FDR),])))
    print(paste("DEGs at FDR ≤", WALD_FDR, "& |LFC| ≥", WALD_lfc_cutoff, ":", nrow(res_LIVER_HEAT[which(res_LIVER_HEAT$padj <= WALD_FDR & abs(res_LIVER_HEAT$log2FoldChange) >= WALD_lfc_cutoff),])))
    summary(res_LIVER_HEAT)
    write.table(res_LIVER_HEAT, "Results/Wald_Results/LIVER_HEAT_vs_CONTROL_WALD_DE.tsv", sep="\t", col.names=NA, quote = FALSE)
  } else {
    print("res_LIVER_HEAT not found or is NULL")
  }
  
  # TESTIS_COLD
  if(exists("res_TESTIS_COLD") && !is.null(res_TESTIS_COLD)) {
    print(paste("DEGs at FDR ≤", WALD_FDR, ":", nrow(res_TESTIS_COLD[which(res_TESTIS_COLD$padj <= WALD_FDR),])))
    print(paste("DEGs at FDR ≤", WALD_FDR, "& |LFC| ≥", WALD_lfc_cutoff, ":", nrow(res_TESTIS_COLD[which(res_TESTIS_COLD$padj <= WALD_FDR & abs(res_TESTIS_COLD$log2FoldChange) >= WALD_lfc_cutoff),])))
    summary(res_TESTIS_COLD)
    write.table(res_TESTIS_COLD, "Results/Wald_Results/TESTIS_COLD_vs_CONTROL_WALD_DE.tsv", sep="\t", col.names=NA, quote = FALSE)
  } else {
    print("res_TESTIS_COLD not found or is NULL")
  }
  
  # TESTIS_HEAT
  if(exists("res_TESTIS_HEAT") && !is.null(res_TESTIS_HEAT)) {
    print(paste("DEGs at FDR ≤", WALD_FDR, ":", nrow(res_TESTIS_HEAT[which(res_TESTIS_HEAT$padj <= WALD_FDR),])))
    print(paste("DEGs at FDR ≤", WALD_FDR, "& |LFC| ≥", WALD_lfc_cutoff, ":", nrow(res_TESTIS_HEAT[which(res_TESTIS_HEAT$padj <= WALD_FDR & abs(res_TESTIS_HEAT$log2FoldChange) >= WALD_lfc_cutoff),])))
    summary(res_TESTIS_HEAT)
    write.table(res_TESTIS_HEAT, "Results/Wald_Results/TESTIS_HEAT_vs_CONTROL_WALD_DE.tsv", sep="\t", col.names=NA, quote = FALSE)
  } else {
    print("res_TESTIS_HEAT not found or is NULL")
  }
}

#######################
# Summary Table
#######################
{
  # Collect DEG counts
  deg_count <- function(res_object) {
    if(exists(deparse(substitute(res_object))) && !is.null(res_object)) {
      return(nrow(res_object[which(res_object$padj <= WALD_FDR),]))
    } else {
      return(0)
    }
  }
  
  deg_lfc_count <- function(res_object) {
    if(exists(deparse(substitute(res_object))) && !is.null(res_object)) {
      return(nrow(res_object[which(res_object$padj <= WALD_FDR & abs(res_object$log2FoldChange) >= WALD_lfc_cutoff),]))
    } else {
      return(0)
    }
  }
  
  # Collect DEG counts in specified order
  wald_deg_counts <- c(
    LIVER_HEAT = deg_count(res_LIVER_HEAT),
    LIVER_COLD = deg_count(res_LIVER_COLD),
    TESTIS_HEAT = deg_count(res_TESTIS_HEAT),
    TESTIS_COLD = deg_count(res_TESTIS_COLD),
    BRAIN_HEAT = deg_count(res_BRAIN_HEAT),
    BRAIN_COLD = deg_count(res_BRAIN_COLD),
    HEART_HEAT = deg_count(res_HEART_HEAT),
    HEART_COLD = deg_count(res_HEART_COLD)
  )
  
  # Collect DEG counts with fold change cutoff in same order
  wald_deg_lfc_counts <- c(
    LIVER_HEAT = deg_lfc_count(res_LIVER_HEAT),
    LIVER_COLD = deg_lfc_count(res_LIVER_COLD),
    TESTIS_HEAT = deg_lfc_count(res_TESTIS_HEAT),
    TESTIS_COLD = deg_lfc_count(res_TESTIS_COLD),
    BRAIN_HEAT = deg_lfc_count(res_BRAIN_HEAT),
    BRAIN_COLD = deg_lfc_count(res_BRAIN_COLD),
    HEART_HEAT = deg_lfc_count(res_HEART_HEAT),
    HEART_COLD = deg_lfc_count(res_HEART_COLD)
  )
  
  # Create summary table
  wald_summary_table <- data.frame(
    TissueComparison = names(wald_deg_counts),
    DEGs_FDR = as.numeric(wald_deg_counts),
    DEGs_FDR_LFC = as.numeric(wald_deg_lfc_counts),
    stringsAsFactors = FALSE
  )
  
  # Add tissue and comparison columns for easier analysis
  wald_summary_table$Tissue <- gsub("_(COLD|HEAT)$", "", wald_summary_table$TissueComparison)
  wald_summary_table$Comparison <- gsub("^.*_", "", wald_summary_table$TissueComparison)
  
  # Don't reorder - keep the specified order
  print("WALD Summary Table (in specified order):")
  print(wald_summary_table)
  
  # Save summary
  write.table(wald_summary_table, "Results/Wald_Results/WALD_DEG_Summary.tsv", 
              quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
}

#######################
# Write L2FC SHRUNKEN Results Tables
#######################
{
  # BRAIN_COLD
  res_WALD_BRAIN_COLD_shrink <- lfcShrink(dds_BRAIN_COLD, contrast = c("CONDITION", "COLD", "CONTROL"), type="ashr", parallel = TRUE)
  res_WALD_BRAIN_COLD_shrink <- res_WALD_BRAIN_COLD_shrink[which(!is.na(res_WALD_BRAIN_COLD_shrink$pvalue)),]
  res_WALD_BRAIN_COLD_shrink <- res_WALD_BRAIN_COLD_shrink[order(res_WALD_BRAIN_COLD_shrink$padj),]
  write.table(res_WALD_BRAIN_COLD_shrink, "Results/Wald_Results/ShrunkenL2FC/BRAIN_COLD_vs_CONTROL_WALD_Shrunk.tsv", sep="\t", col.names=NA, quote=FALSE)
  
  # BRAIN_HEAT
  res_WALD_BRAIN_HEAT_shrink <- lfcShrink(dds_BRAIN_HEAT, contrast = c("CONDITION", "HEAT", "CONTROL"), type="ashr", parallel = TRUE)
  res_WALD_BRAIN_HEAT_shrink <- res_WALD_BRAIN_HEAT_shrink[which(!is.na(res_WALD_BRAIN_HEAT_shrink$pvalue)),]
  res_WALD_BRAIN_HEAT_shrink <- res_WALD_BRAIN_HEAT_shrink[order(res_WALD_BRAIN_HEAT_shrink$padj),]
  write.table(res_WALD_BRAIN_HEAT_shrink, "Results/Wald_Results/ShrunkenL2FC/BRAIN_HEAT_vs_CONTROL_WALD_Shrunk.tsv", sep="\t", col.names=NA, quote=FALSE)
  
  # HEART_COLD
  res_WALD_HEART_COLD_shrink <- lfcShrink(dds_HEART_COLD, contrast = c("CONDITION", "COLD", "CONTROL"), type="ashr", parallel = TRUE)
  res_WALD_HEART_COLD_shrink <- res_WALD_HEART_COLD_shrink[which(!is.na(res_WALD_HEART_COLD_shrink$pvalue)),]
  res_WALD_HEART_COLD_shrink <- res_WALD_HEART_COLD_shrink[order(res_WALD_HEART_COLD_shrink$padj),]
  write.table(res_WALD_HEART_COLD_shrink, "Results/Wald_Results/ShrunkenL2FC/HEART_COLD_vs_CONTROL_WALD_Shrunk.tsv", sep="\t", col.names=NA, quote=FALSE)
  
  # HEART_HEAT
  res_WALD_HEART_HEAT_shrink <- lfcShrink(dds_HEART_HEAT, contrast = c("CONDITION", "HEAT", "CONTROL"), type="ashr", parallel = TRUE)
  res_WALD_HEART_HEAT_shrink <- res_WALD_HEART_HEAT_shrink[which(!is.na(res_WALD_HEART_HEAT_shrink$pvalue)),]
  res_WALD_HEART_HEAT_shrink <- res_WALD_HEART_HEAT_shrink[order(res_WALD_HEART_HEAT_shrink$padj),]
  write.table(res_WALD_HEART_HEAT_shrink, "Results/Wald_Results/ShrunkenL2FC/HEART_HEAT_vs_CONTROL_WALD_Shrunk.tsv", sep="\t", col.names=NA, quote=FALSE)
  
  # LIVER_COLD
  res_WALD_LIVER_COLD_shrink <- lfcShrink(dds_LIVER_COLD, contrast = c("CONDITION", "COLD", "CONTROL"), type="ashr", parallel = TRUE)
  res_WALD_LIVER_COLD_shrink <- res_WALD_LIVER_COLD_shrink[which(!is.na(res_WALD_LIVER_COLD_shrink$pvalue)),]
  res_WALD_LIVER_COLD_shrink <- res_WALD_LIVER_COLD_shrink[order(res_WALD_LIVER_COLD_shrink$padj),]
  write.table(res_WALD_LIVER_COLD_shrink, "Results/Wald_Results/ShrunkenL2FC/LIVER_COLD_vs_CONTROL_WALD_Shrunk.tsv", sep="\t", col.names=NA, quote=FALSE)
  
  # LIVER_HEAT
  res_WALD_LIVER_HEAT_shrink <- lfcShrink(dds_LIVER_HEAT, contrast = c("CONDITION", "HEAT", "CONTROL"), type="ashr", parallel = TRUE)
  res_WALD_LIVER_HEAT_shrink <- res_WALD_LIVER_HEAT_shrink[which(!is.na(res_WALD_LIVER_HEAT_shrink$pvalue)),]
  res_WALD_LIVER_HEAT_shrink <- res_WALD_LIVER_HEAT_shrink[order(res_WALD_LIVER_HEAT_shrink$padj),]
  write.table(res_WALD_LIVER_HEAT_shrink, "Results/Wald_Results/ShrunkenL2FC/LIVER_HEAT_vs_CONTROL_WALD_Shrunk.tsv", sep="\t", col.names=NA, quote=FALSE)
  
  # TESTIS_COLD
  res_WALD_TESTIS_COLD_shrink <- lfcShrink(dds_TESTIS_COLD, contrast = c("CONDITION", "COLD", "CONTROL"), type="ashr", parallel = TRUE)
  res_WALD_TESTIS_COLD_shrink <- res_WALD_TESTIS_COLD_shrink[which(!is.na(res_WALD_TESTIS_COLD_shrink$pvalue)),]
  res_WALD_TESTIS_COLD_shrink <- res_WALD_TESTIS_COLD_shrink[order(res_WALD_TESTIS_COLD_shrink$padj),]
  write.table(res_WALD_TESTIS_COLD_shrink, "Results/Wald_Results/ShrunkenL2FC/TESTIS_COLD_vs_CONTROL_WALD_Shrunk.tsv", sep="\t", col.names=NA, quote=FALSE)
  
  # TESTIS_HEAT
  res_WALD_TESTIS_HEAT_shrink <- lfcShrink(dds_TESTIS_HEAT, contrast = c("CONDITION", "HEAT", "CONTROL"), type="ashr", parallel = TRUE)
  res_WALD_TESTIS_HEAT_shrink <- res_WALD_TESTIS_HEAT_shrink[which(!is.na(res_WALD_TESTIS_HEAT_shrink$pvalue)),]
  res_WALD_TESTIS_HEAT_shrink <- res_WALD_TESTIS_HEAT_shrink[order(res_WALD_TESTIS_HEAT_shrink$padj),]
  write.table(res_WALD_TESTIS_HEAT_shrink, "Results/Wald_Results/ShrunkenL2FC/TESTIS_HEAT_vs_CONTROL_WALD_Shrunk.tsv", sep="\t", col.names=NA, quote=FALSE)
}


#######################
# Extract Normalized Counts
#######################
{
  # Extract normalized counts for each pairwise comparison
  norm_counts_BRAIN_COLD <- counts(dds_BRAIN_COLD, normalized=TRUE)
  norm_counts_BRAIN_HEAT <- counts(dds_BRAIN_HEAT, normalized=TRUE)
  norm_counts_HEART_COLD <- counts(dds_HEART_COLD, normalized=TRUE)
  norm_counts_HEART_HEAT <- counts(dds_HEART_HEAT, normalized=TRUE)
  norm_counts_LIVER_COLD <- counts(dds_LIVER_COLD, normalized=TRUE)
  norm_counts_LIVER_HEAT <- counts(dds_LIVER_HEAT, normalized=TRUE)
  norm_counts_TESTIS_COLD <- counts(dds_TESTIS_COLD, normalized=TRUE)
  norm_counts_TESTIS_HEAT <- counts(dds_TESTIS_HEAT, normalized=TRUE)
  
  # Save normalized counts for each comparison
  write.table(norm_counts_BRAIN_COLD, "Results/Wald_Results/Normalized_Counts/BRAIN_COLD_normalized_counts.tsv", sep="\t", col.names=NA, quote=FALSE)
  write.table(norm_counts_BRAIN_HEAT, "Results/Wald_Results/Normalized_Counts/BRAIN_HEAT_normalized_counts.tsv", sep="\t", col.names=NA, quote=FALSE)
  write.table(norm_counts_HEART_COLD, "Results/Wald_Results/Normalized_Counts/HEART_COLD_normalized_counts.tsv", sep="\t", col.names=NA, quote=FALSE)
  write.table(norm_counts_HEART_HEAT, "Results/Wald_Results/Normalized_Counts/HEART_HEAT_normalized_counts.tsv", sep="\t", col.names=NA, quote=FALSE)
  write.table(norm_counts_LIVER_COLD, "Results/Wald_Results/Normalized_Counts/LIVER_COLD_normalized_counts.tsv", sep="\t", col.names=NA, quote=FALSE)
  write.table(norm_counts_LIVER_HEAT, "Results/Wald_Results/Normalized_Counts/LIVER_HEAT_normalized_counts.tsv", sep="\t", col.names=NA, quote=FALSE)
  write.table(norm_counts_TESTIS_COLD, "Results/Wald_Results/Normalized_Counts/TESTIS_COLD_normalized_counts.tsv", sep="\t", col.names=NA, quote=FALSE)
  write.table(norm_counts_TESTIS_HEAT, "Results/Wald_Results/Normalized_Counts/TESTIS_HEAT_normalized_counts.tsv", sep="\t", col.names=NA, quote=FALSE)
  
  print("Normalized counts extracted and saved for all pairwise comparisons")
}

#######################
# WALD Heatmaps
#######################

# BRAIN_COLD
DE_WALD_BRAIN_COLD <- rownames(res_BRAIN_COLD[which(res_BRAIN_COLD$padj <= WALD_FDR),])  
print(paste("BRAIN_COLD DE genes for heatmap:", length(DE_WALD_BRAIN_COLD)))
if(length(DE_WALD_BRAIN_COLD) > 0) {
  DE_WALD_BRAIN_COLD_matrix <- counts(dds_BRAIN_COLD, normalized=FALSE)[DE_WALD_BRAIN_COLD,]
  vst_DE_WALD_BRAIN_COLD <- varianceStabilizingTransformation(DE_WALD_BRAIN_COLD_matrix, blind = TRUE, fitType = "parametric")
  png(filename="Heatmaps/Wald_Heatmaps/BRAIN_COLD_vs_CONTROL_WALD_heatmap.png", units="in", width=12, height=13, pointsize=12, res=800)
  heatmap3(vst_DE_WALD_BRAIN_COLD, method = "complete", Rowv=TRUE, Colv=NA, col=hmcol, scale="row", labRow=NA, showRowDendro = TRUE)
  dev.off()
} else {
  print("No DE genes for BRAIN_COLD heatmap")
}

# BRAIN_HEAT
DE_WALD_BRAIN_HEAT <- rownames(res_BRAIN_HEAT[which(res_BRAIN_HEAT$padj <= WALD_FDR),])  
print(paste("BRAIN_HEAT DE genes for heatmap:", length(DE_WALD_BRAIN_HEAT)))
if(length(DE_WALD_BRAIN_HEAT) > 0) {
  DE_WALD_BRAIN_HEAT_matrix <- counts(dds_BRAIN_HEAT, normalized=FALSE)[DE_WALD_BRAIN_HEAT,]
  vst_DE_WALD_BRAIN_HEAT <- varianceStabilizingTransformation(DE_WALD_BRAIN_HEAT_matrix, blind = TRUE, fitType = "parametric")
  png(filename="Heatmaps/Wald_Heatmaps/BRAIN_HEAT_vs_CONTROL_WALD_heatmap.png", units="in", width=12, height=13, pointsize=12, res=800)
  heatmap3(vst_DE_WALD_BRAIN_HEAT, method = "complete", Rowv=TRUE, Colv=NA, col=hmcol, scale="row", labRow=NA, showRowDendro = TRUE)
  dev.off()
} else {
  print("No DE genes for BRAIN_HEAT heatmap")
}

# HEART_COLD
DE_WALD_HEART_COLD <- rownames(res_HEART_COLD[which(res_HEART_COLD$padj <= WALD_FDR),])  
print(paste("HEART_COLD DE genes for heatmap:", length(DE_WALD_HEART_COLD)))
if(length(DE_WALD_HEART_COLD) > 0) {
  DE_WALD_HEART_COLD_matrix <- counts(dds_HEART_COLD, normalized=FALSE)[DE_WALD_HEART_COLD,]
  vst_DE_WALD_HEART_COLD <- varianceStabilizingTransformation(DE_WALD_HEART_COLD_matrix, blind = TRUE, fitType = "parametric")
  png(filename="Heatmaps/Wald_Heatmaps/HEART_COLD_vs_CONTROL_WALD_heatmap.png", units="in", width=12, height=13, pointsize=12, res=800)
  heatmap3(vst_DE_WALD_HEART_COLD, method = "complete", Rowv=TRUE, Colv=NA, col=hmcol, scale="row", labRow=NA, showRowDendro = TRUE)
  dev.off()
} else {
  print("No DE genes for HEART_COLD heatmap")
}

# HEART_HEAT
DE_WALD_HEART_HEAT <- rownames(res_HEART_HEAT[which(res_HEART_HEAT$padj <= WALD_FDR),])  
print(paste("HEART_HEAT DE genes for heatmap:", length(DE_WALD_HEART_HEAT)))
if(length(DE_WALD_HEART_HEAT) > 0) {
  DE_WALD_HEART_HEAT_matrix <- counts(dds_HEART_HEAT, normalized=FALSE)[DE_WALD_HEART_HEAT,]
  vst_DE_WALD_HEART_HEAT <- varianceStabilizingTransformation(DE_WALD_HEART_HEAT_matrix, blind = TRUE, fitType = "parametric")
  png(filename="Heatmaps/Wald_Heatmaps/HEART_HEAT_vs_CONTROL_WALD_heatmap.png", units="in", width=12, height=13, pointsize=12, res=800)
  heatmap3(vst_DE_WALD_HEART_HEAT, method = "complete", Rowv=TRUE, Colv=NA, col=hmcol, scale="row", labRow=NA, showRowDendro = TRUE)
  dev.off()
} else {
  print("No DE genes for HEART_HEAT heatmap")
}

# LIVER_COLD
DE_WALD_LIVER_COLD <- rownames(res_LIVER_COLD[which(res_LIVER_COLD$padj <= WALD_FDR),])  
print(paste("LIVER_COLD DE genes for heatmap:", length(DE_WALD_LIVER_COLD)))
if(length(DE_WALD_LIVER_COLD) > 0) {
  DE_WALD_LIVER_COLD_matrix <- counts(dds_LIVER_COLD, normalized=FALSE)[DE_WALD_LIVER_COLD,]
  vst_DE_WALD_LIVER_COLD <- varianceStabilizingTransformation(DE_WALD_LIVER_COLD_matrix, blind = TRUE, fitType = "parametric")
  png(filename="Heatmaps/Wald_Heatmaps/LIVER_COLD_vs_CONTROL_WALD_heatmap.png", units="in", width=12, height=13, pointsize=12, res=800)
  heatmap3(vst_DE_WALD_LIVER_COLD, method = "complete", Rowv=TRUE, Colv=NA, col=hmcol, scale="row", labRow=NA, showRowDendro = TRUE)
  dev.off()
} else {
  print("No DE genes for LIVER_COLD heatmap")
}

# LIVER_HEAT
DE_WALD_LIVER_HEAT <- rownames(res_LIVER_HEAT[which(res_LIVER_HEAT$padj <= WALD_FDR),])  
print(paste("LIVER_HEAT DE genes for heatmap:", length(DE_WALD_LIVER_HEAT)))
if(length(DE_WALD_LIVER_HEAT) > 0) {
  DE_WALD_LIVER_HEAT_matrix <- counts(dds_LIVER_HEAT, normalized=FALSE)[DE_WALD_LIVER_HEAT,]
  vst_DE_WALD_LIVER_HEAT <- varianceStabilizingTransformation(DE_WALD_LIVER_HEAT_matrix, blind = TRUE, fitType = "parametric")
  png(filename="Heatmaps/Wald_Heatmaps/LIVER_HEAT_vs_CONTROL_WALD_heatmap.png", units="in", width=12, height=13, pointsize=12, res=800)
  heatmap3(vst_DE_WALD_LIVER_HEAT, method = "complete", Rowv=TRUE, Colv=NA, col=hmcol, scale="row", labRow=NA, showRowDendro = TRUE)
  dev.off()
} else {
  print("No DE genes for LIVER_HEAT heatmap")
}

# TESTIS_COLD
DE_WALD_TESTIS_COLD <- rownames(res_TESTIS_COLD[which(res_TESTIS_COLD$padj <= WALD_FDR),])  
print(paste("TESTIS_COLD DE genes for heatmap:", length(DE_WALD_TESTIS_COLD)))
if(length(DE_WALD_TESTIS_COLD) > 0) {
  DE_WALD_TESTIS_COLD_matrix <- counts(dds_TESTIS_COLD, normalized=FALSE)[DE_WALD_TESTIS_COLD,]
  vst_DE_WALD_TESTIS_COLD <- varianceStabilizingTransformation(DE_WALD_TESTIS_COLD_matrix, blind = TRUE, fitType = "parametric")
  png(filename="Heatmaps/Wald_Heatmaps/TESTIS_COLD_vs_CONTROL_WALD_heatmap.png", units="in", width=12, height=13, pointsize=12, res=800)
  heatmap3(vst_DE_WALD_TESTIS_COLD, method = "complete", Rowv=TRUE, Colv=NA, col=hmcol, scale="row", labRow=NA, showRowDendro = TRUE)
  dev.off()
} else {
  print("No DE genes for TESTIS_COLD heatmap")
}

# TESTIS_HEAT
DE_WALD_TESTIS_HEAT <- rownames(res_TESTIS_HEAT[which(res_TESTIS_HEAT$padj <= WALD_FDR),])  
print(paste("TESTIS_HEAT DE genes for heatmap:", length(DE_WALD_TESTIS_HEAT)))
if(length(DE_WALD_TESTIS_HEAT) > 0) {
  DE_WALD_TESTIS_HEAT_matrix <- counts(dds_TESTIS_HEAT, normalized=FALSE)[DE_WALD_TESTIS_HEAT,]
  vst_DE_WALD_TESTIS_HEAT <- varianceStabilizingTransformation(DE_WALD_TESTIS_HEAT_matrix, blind = TRUE, fitType = "parametric")
  png(filename="Heatmaps/Wald_Heatmaps/TESTIS_HEAT_vs_CONTROL_WALD_heatmap.png", units="in", width=12, height=13, pointsize=12, res=800)
  heatmap3(vst_DE_WALD_TESTIS_HEAT, method = "complete", Rowv=TRUE, Colv=NA, col=hmcol, scale="row", labRow=NA, showRowDendro = TRUE)
  dev.off()
} else {
  print("No DE genes for TESTIS_HEAT heatmap")
}


#######################
# WALD Volcano plots
#######################

# BRAIN_COLD
png(filename="Volcano_Plots/BRAIN_COLD_vs_CONTROL_WALD_volcano.png", units="in", width=10, height=10, pointsize=4, res=800)
EnhancedVolcano(res_BRAIN_COLD,
                lab = rownames(res_BRAIN_COLD), 
                x = 'log2FoldChange',
                y = 'padj',
                #xlim = c(-10, 10),
                pointSize = 4,
                labSize = 0,
                pCutoff = 0.05,
                FCcutoff=WALD_lfc_cutoff,
                gridlines.major = TRUE,
                gridlines.minor = FALSE,
                col=c('black', 'blue', 'black', 'darkorange'),
                title = 'Differential expression in Brain: Cold vs Control',
                titleLabSize = 14,
                subtitle ="",
                legendPosition = 'none',
                caption = "Orange: FDR <= 0.05 & raw fold change > 1.5x",
                border = 'full',
                borderWidth = 0.5,
                borderColour = 'black',
)
dev.off()

# BRAIN_HEAT
png(filename="Volcano_Plots/BRAIN_HEAT_vs_CONTROL_WALD_volcano.png", units="in", width=10, height=10, pointsize=4, res=800)
EnhancedVolcano(res_BRAIN_HEAT,
                lab = rownames(res_BRAIN_HEAT), 
                x = 'log2FoldChange',
                y = 'padj',
                #xlim = c(-10, 10),
                pointSize = 4,
                labSize = 0,
                pCutoff = 0.05,
                FCcutoff=WALD_lfc_cutoff,
                gridlines.major = TRUE,
                gridlines.minor = FALSE,
                col=c('black', 'blue', 'black', 'darkorange'),
                title = 'Differential expression in Brain: Heat vs Control',
                titleLabSize = 14,
                subtitle ="",
                legendPosition = 'none',
                caption = "Orange: FDR <= 0.05 & raw fold change > 1.5x",
                border = 'full',
                borderWidth = 0.5,
                borderColour = 'black',
)
dev.off()

# HEART_COLD
png(filename="Volcano_Plots/HEART_COLD_vs_CONTROL_WALD_volcano.png", units="in", width=10, height=10, pointsize=4, res=800)
EnhancedVolcano(res_HEART_COLD,
                lab = rownames(res_HEART_COLD), 
                x = 'log2FoldChange',
                y = 'padj',
                #xlim = c(-10, 10),
                pointSize = 4,
                labSize = 0,
                pCutoff = 0.05,
                FCcutoff=WALD_lfc_cutoff,
                gridlines.major = TRUE,
                gridlines.minor = FALSE,
                col=c('black', 'blue', 'black', 'darkorange'),
                title = 'Differential expression in Heart: Cold vs Control',
                titleLabSize = 14,
                subtitle ="",
                legendPosition = 'none',
                caption = "Orange: FDR <= 0.05 & raw fold change > 1.5x",
                border = 'full',
                borderWidth = 0.5,
                borderColour = 'black',
)
dev.off()

# HEART_HEAT
png(filename="Volcano_Plots/HEART_HEAT_vs_CONTROL_WALD_volcano.png", units="in", width=10, height=10, pointsize=4, res=800)
EnhancedVolcano(res_HEART_HEAT,
                lab = rownames(res_HEART_HEAT), 
                x = 'log2FoldChange',
                y = 'padj',
                #xlim = c(-10, 10),
                pointSize = 4,
                labSize = 0,
                pCutoff = 0.05,
                FCcutoff=WALD_lfc_cutoff,
                gridlines.major = TRUE,
                gridlines.minor = FALSE,
                col=c('black', 'blue', 'black', 'darkorange'),
                title = 'Differential expression in Heart: Heat vs Control',
                titleLabSize = 14,
                subtitle ="",
                legendPosition = 'none',
                caption = "Orange: FDR <= 0.05 & raw fold change > 1.5x",
                border = 'full',
                borderWidth = 0.5,
                borderColour = 'black',
)
dev.off()

# LIVER_COLD
png(filename="Volcano_Plots/LIVER_COLD_vs_CONTROL_WALD_volcano.png", units="in", width=10, height=10, pointsize=4, res=800)
EnhancedVolcano(res_LIVER_COLD,
                lab = rownames(res_LIVER_COLD), 
                x = 'log2FoldChange',
                y = 'padj',
                #xlim = c(-10, 10),
                pointSize = 4,
                labSize = 0,
                pCutoff = 0.05,
                FCcutoff=WALD_lfc_cutoff,
                gridlines.major = TRUE,
                gridlines.minor = FALSE,
                col=c('black', 'blue', 'black', 'darkorange'),
                title = 'Differential expression in Liver: Cold vs Control',
                titleLabSize = 14,
                subtitle ="",
                legendPosition = 'none',
                caption = "Orange: FDR <= 0.05 & raw fold change > 1.5x",
                border = 'full',
                borderWidth = 0.5,
                borderColour = 'black',
)
dev.off()

# LIVER_HEAT
png(filename="Volcano_Plots/LIVER_HEAT_vs_CONTROL_WALD_volcano.png", units="in", width=10, height=10, pointsize=4, res=800)
EnhancedVolcano(res_LIVER_HEAT,
                lab = rownames(res_LIVER_HEAT), 
                x = 'log2FoldChange',
                y = 'padj',
                #xlim = c(-10, 10),
                pointSize = 4,
                labSize = 0,
                pCutoff = 0.05,
                FCcutoff=WALD_lfc_cutoff,
                gridlines.major = TRUE,
                gridlines.minor = FALSE,
                col=c('black', 'blue', 'black', 'darkorange'),
                title = 'Differential expression in Liver: Heat vs Control',
                titleLabSize = 14,
                subtitle ="",
                legendPosition = 'none',
                caption = "Orange: FDR <= 0.05 & raw fold change > 1.5x",
                border = 'full',
                borderWidth = 0.5,
                borderColour = 'black',
)
dev.off()

# TESTIS_COLD
png(filename="Volcano_Plots/TESTIS_COLD_vs_CONTROL_WALD_volcano.png", units="in", width=10, height=10, pointsize=4, res=800)
EnhancedVolcano(res_TESTIS_COLD,
                lab = rownames(res_TESTIS_COLD), 
                x = 'log2FoldChange',
                y = 'padj',
                #xlim = c(-10, 10),
                pointSize = 4,
                labSize = 0,
                pCutoff = 0.05,
                FCcutoff=WALD_lfc_cutoff,
                gridlines.major = TRUE,
                gridlines.minor = FALSE,
                col=c('black', 'blue', 'black', 'darkorange'),
                title = 'Differential expression in Testis: Cold vs Control',
                titleLabSize = 14,
                subtitle ="",
                legendPosition = 'none',
                caption = "Orange: FDR <= 0.05 & raw fold change > 1.5x",
                border = 'full',
                borderWidth = 0.5,
                borderColour = 'black',
)
dev.off()

# TESTIS_HEAT
png(filename="Volcano_Plots/TESTIS_HEAT_vs_CONTROL_WALD_volcano.png", units="in", width=10, height=10, pointsize=4, res=800)
EnhancedVolcano(res_TESTIS_HEAT,
                lab = rownames(res_TESTIS_HEAT), 
                x = 'log2FoldChange',
                y = 'padj',
                #xlim = c(-10, 10),
                pointSize = 4,
                labSize = 0,
                pCutoff = 0.05,
                FCcutoff=WALD_lfc_cutoff,
                gridlines.major = TRUE,
                gridlines.minor = FALSE,
                col=c('black', 'blue', 'black', 'darkorange'),
                title = 'Differential expression in Testis: Heat vs Control',
                titleLabSize = 14,
                subtitle ="",
                legendPosition = 'none',
                caption = "Orange: FDR <= 0.05 & raw fold change > 1.5x",
                border = 'full',
                borderWidth = 0.5,
                borderColour = 'black',
)
dev.off()

