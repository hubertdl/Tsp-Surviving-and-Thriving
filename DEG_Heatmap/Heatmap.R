# the following will generate a heat of the top differentially expressed genes by tissue

setwd("")

library(tidyverse)

# Clear environment
rm(list = ls())+
graphics.off()+
cat("\014")

# Load and filter your data as before
df <- read.csv("DE_ALL.csv")
head(df)

ref <- read.csv("Tsp_transcript_info_updated.csv")
head(ref)
str(ref)

df_significant <- df %>%
  group_by(tissue, ID) %>%
  filter(any(!is.na(padj))) %>%  # keep only groups with at least one non-NA padj
  summarise(min_padj = min(padj, na.rm = TRUE), .groups = "drop") %>%
  group_by(tissue) %>%
  slice_min(order_by = min_padj, n = 20, with_ties = FALSE) %>%
  left_join(df, by = c("tissue", "ID"))
nrow(df_significant)
head(df_significant)

# Pivot to wide format for clustering
df_wide <- df_significant %>%
  select(tissue, treatment, ID, L2FC) %>%
  pivot_wider(names_from = treatment, values_from = L2FC)
head(df_wide)

# Order by expression level within each tissue (highest to lowest)
clustered_df <- df_wide %>%
  group_by(tissue) %>%
  mutate(
    # Calculate mean expression across treatments for ordering
    mean_expression = rowMeans(cbind(cold, hot), na.rm = TRUE),
    # Create ordering (highest expression first)
    expression_order = rank(-mean_expression, ties.method = "first")
  ) %>%
  arrange(tissue, expression_order) %>%
  ungroup()

# Merge back into long format for ggplot
df_plot <- clustered_df %>%
  pivot_longer(cols = -c(tissue, ID, mean_expression, expression_order), names_to = "treatment", values_to = "L2FC")

# Add gene names
df_plot$gene <- ref$gene[match(df_plot$ID, ref$ID)]
df_plot$gene <- sub("_$", "", df_plot$gene)
df_plot$gene <- gsub("_", " ", df_plot$gene)
df_plot$gene <- ifelse(df_plot$gene == "", as.character(df_plot$ID), df_plot$gene)

# Plot
#png("Transcript_LFC_by_treatment_and_Tissue.png", width = 12, height = 8, units = "in", res = 1000)
# Plot with expression-based ordering
ggplot(df_plot, aes(x = treatment, y = reorder(ID, -expression_order), fill = L2FC)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#007aa5", mid = "white", high = "#b22222", midpoint = 0,
    limits = c(-5, 5), oob = scales::squish, name = "L2FC"
  ) +
  scale_y_discrete(labels = setNames(df_plot$gene, df_plot$ID)) +
  facet_wrap(~ tissue, ncol = 2, scales = "free_y") +
  labs(title = "",
       x = "Treatment", y = "Transcript") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8, hjust = 1),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1, size = 12),
    strip.text = element_text(size = 14)
  )

#dev.off()


# Extract unique genes from the plot data
genes_in_plot <- df_plot %>%
  select(ID, gene) %>%
  distinct() %>%
  arrange(ID)

# Write to CSV file
write.csv(genes_in_plot, "genes_in_heatmap.csv", row.names = FALSE)

