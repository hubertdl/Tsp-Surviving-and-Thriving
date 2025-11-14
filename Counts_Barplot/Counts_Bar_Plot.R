setwd("")

library(ggplot2)
library(reshape2)
library(forcats)
library(gridExtra)
library(grid)

# Clear environment, graphics, and console
rm(list = ls())+
graphics.off()+
cat("/014")

# Read in counts of DE genes
df <- read.table("DE_Counts.tsv", header = TRUE, sep = "\t")
df <- melt(df)
df <- as.data.frame(df)

# Create separate grouping variables
df$TREATMENT <- ifelse(grepl("^Cold", df$variable), "COLD", "HEAT")
df$DIRECTION <- ifelse(grepl("Up$", df$variable), "UP", "DOWN")
df$TISSUE <- df$Tissue
df$GROUP <- paste(df$TREATMENT, df$TISSUE)
df$TREATMENT_DIRECTION <- paste(df$TREATMENT, df$DIRECTION, sep = "_")

# Create capped values for plotting - cap at chart limits
df$value_capped <- ifelse(df$value > 450, 450, ifelse(df$value < -450, -450, df$value))
df$is_capped <- abs(df$value) > 450
# Split data by treatment
heat_df <- df[df$TREATMENT == "HEAT", ]
cold_df <- df[df$TREATMENT == "COLD", ]

# Set tissue order
tissue_order <- rev(c("Brain", "Heart", "Liver", "Testis"))
heat_df$TISSUE <- factor(heat_df$TISSUE, levels = tissue_order)
cold_df$TISSUE <- factor(cold_df$TISSUE, levels = tissue_order)

# Define fill colors
heat_colors <- scale_fill_manual(values = c("HEAT_UP" = "#b22222", "HEAT_DOWN" = "#007aa5"))
cold_colors <- scale_fill_manual(values = c("COLD_UP" = "#b22222", "COLD_DOWN" = "#007aa5"))

# Heat plot
heat_plot <- ggplot(heat_df, aes(y = TISSUE, x = value_capped, fill = TREATMENT_DIRECTION)) +
  geom_col(colour = "black", width = .95) +
  coord_cartesian(xlim = c(-55, 70)) +
  scale_x_continuous("", expand = c(0, 0)) +
  scale_y_discrete("") +
  heat_colors +
  geom_text(aes(label = ifelse(value < 0, -value, value)), 
            hjust = ifelse(heat_df$value < 0, 1.2, ifelse(heat_df$value == 0, 1.2, -0.2))) +
  annotate("text", x = -40, y = 4.5, label = "Decreased", angle = 0, vjust = 5.0, hjust = 0.5, size = 5) +
  annotate("text", x = 60, y = 4.5, label = "Increased", angle = 0, vjust = 5.0, hjust = 2.0, size = 5) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y = element_text(hjust = 1, size = 14),
    axis.text.x = element_blank(),
    axis.title.x = element_text(hjust = 0.425),
    legend.title = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.425)
  ) +
  ggtitle("Heat Stress")

# Cold plot
cold_plot <- ggplot(cold_df, aes(y = TISSUE, x = value_capped, fill = TREATMENT_DIRECTION)) +
  geom_col(colour = "black", width = .95) +
  coord_cartesian(xlim = c(-55, 70)) +
  scale_x_continuous("Number of D.E. Genes", expand = c(0, 0)) +
  scale_y_discrete("") +
  cold_colors +
  geom_text(aes(label = ifelse(value < 0, -value, value)), 
            hjust = ifelse(cold_df$value < 0, 1.2, ifelse(cold_df$value == 0, 1.2, -0.2))) +
  annotate("text", x = -40, y = 4.5, label = "Decreased", angle = 0, vjust = 5.0, hjust = 0.5, size = 5) +
  annotate("text", x = 60, y = 4.5, label = "Increased", angle = 0, vjust = 5.0, hjust = 2.0, size = 5) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y = element_text(hjust = 1, size = 14),
    axis.text.x = element_blank(),
    axis.title.x = element_text(hjust = 0.425),
    legend.title = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.425)
  ) +
  ggtitle("Cold Stress")

# Print combined plots
grid.arrange(heat_plot, cold_plot, nrow = 2)

# Save
png(filename="DE_Counts_barplot.png", units="in", width=10, height=6, res=1000, bg="transparent")
grid.arrange(heat_plot, cold_plot, nrow = 2)
dev.off()
