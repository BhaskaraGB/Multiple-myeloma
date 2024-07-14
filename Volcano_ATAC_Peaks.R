# Clear the workspace
rm(list = ls())

# Load necessary libraries
library(ggplot2)
library(ggthemes)
library(ggrepel)
library(tidyverse)
library(patchwork)

# Set working directory
setwd("/Users/bgovinal/Documents/MDACC/Labmembers/Candy/ATAC_Seq_Ada2b")

# Read input CSV file
atac_data <- read.csv("NT-3d-Dox_vs_s2B-3d-Dox_edgeR_All.-50000-TSS-50000.closest_gene.only_original_peaks_kept.csv",
                      na.strings = c("", "NA"))

# View column names and first few rows
names(atac_data)
head(atac_data)

# Process data: Select relevant columns and categorize differential expression
processed_data <- atac_data %>% 
  select(Gene, log2ratio, FDR) %>% 
  mutate(DiffExp = ifelse(FDR <= 0.05 & abs(log2ratio) >= 0, 
                          ifelse(log2ratio >= 0, "Increased", "Decreased"),
                          "Non-significant"))

# Plot data with ggplot2
volcano_plot <- ggplot(processed_data, aes(log2ratio, -log10(FDR), color = DiffExp)) +
  geom_point() +
  scale_color_manual(values = c("blue", "red", "darkgrey"), name = NULL) +
  geom_vline(xintercept = c(-0.58, 0.58), col = "black", lwd = 0.5, lty = 2) +
  geom_hline(yintercept = -log10(0.05), col = "black", lwd = 0.5, lty = 2) +             
  theme_par() +
  labs(title = "2B vs NT 3d-Dox",
       x = expression ("Log"[2] ~ "FC"),
       y = expression ("-Log"[10] ~ "FDR")) +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        title = element_text(size = 20))

# Display the plot
print(volcano_plot)



# Read input CSV file for 5-day data
atac_data_5day <- read.csv("NT-5d-Dox_vs_s2B-5d-Dox_edgeR_All.-50000-TSS-50000.closest_gene.only_original_peaks_kept.csv", na.strings = c("", "NA"))

# Process data: Select relevant columns and categorize differential expression
processed_data_5day <- atac_data_5day %>% 
  select(Gene, log2ratio, FDR) %>% 
  mutate(DiffExp = ifelse(FDR <= 0.05 & abs(log2ratio) >= 0, 
                          ifelse(log2ratio >= 0, "Increased", "Decreased"),
                          "Non-significant"))

table(processed_data_5day$DiffExp)

# Plot data with ggplot2 for 5-day data
volcano_plot_5day <- ggplot(processed_data_5day, aes(log2ratio, -log10(FDR), color = DiffExp, label = delabel)) +
  geom_point() +
  scale_color_manual(values = c("blue", "red", "darkgrey"), name = NULL) +
  geom_vline(xintercept = c(-0.58, 0.58), col = "black", lwd = 0.5, lty = 2) +
  geom_hline(yintercept = -log10(0.05), col = "black", lwd = 0.5, lty = 2) +             
  theme_par() +
  labs(title = "2B vs NT 5d-Dox",
       x = expression ("Log"[2] ~ "FC"),
       y = expression ("-Log"[10] ~ "FDR")) +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        title = element_text(size = 20),
        legend.text = element_text(size = 20)) +
  guides(color = guide_legend(override.aes = list(size = 8)))

# Display the plot for 5-day data
print(volcano_plot_5day)

# Combine with 3-day plot
volcano_ATAC <- volcano_plot_3day + volcano_plot_5day

# Display the combined plot
print(volcano_ATAC)

# Save the combined plot
ggsave(volcano_ATAC, filename = "Volcano_ATAC.tiff_nonames_Feb7th2024.tiff", dpi = 300, width = 16, height = 8)

