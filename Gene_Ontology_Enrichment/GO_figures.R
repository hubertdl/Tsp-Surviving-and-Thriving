# the following will plot the top significant go terms for each tissue and treatment

setwd("")

#clear environment
rm(list=ls())+
  dev.off()

# read in GO results files
GO <- read.csv("GO_Results_All.csv", header = TRUE)

library(ggplot2)
library(dplyr)
library(tidytext)

# filter for sig and top 20
All_sig <- GO %>%
  filter(!is.na(p.adj), !is.na(delta.rank), p.adj < 0.05) %>%   # sig + clean
  group_by(tissue, condition) %>%
  slice_max(order_by = abs(delta.rank), n = 10, with_ties = FALSE) %>%
  ungroup()

# Global scales for x and bubble size (so visual mappings are comparable)
max_abs   <- max(abs(All_sig$delta.rank), na.rm = TRUE)
max_size  <- max(-log10(All_sig$p.adj), na.rm = TRUE)

# 3) Single faceted plot: rows = tissue, cols = condition → panels align by construction
p_all <- ggplot(
  All_sig,
  aes(x = delta.rank,
      y = reorder_within(name, delta.rank, interaction(tissue, condition)),
      size = -log10(p.adj),
      color = direction)
) +
  geom_point() +
  scale_y_reordered() +
  scale_color_manual(values = c(up = "#b22222", down = "#007aa5")) +
  scale_size(limits = c(0, max_size), range = c(1.5, 6)) +
  coord_cartesian(xlim = c(-max_abs, max_abs)) +
  labs(x = "Delta Rank", y = "GO Term", size = "-log10(FDR)",
       title = "") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    panel.spacing = unit(0.8, "lines")
  ) +
  facet_grid(rows = vars(tissue), cols = vars(condition), scales = "free_y")

# Print or save
p_all


