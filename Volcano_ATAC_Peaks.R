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

# # Label genes based on differential expression
# processed_data$delabel <- NA
# processed_data$delabel[processed_data$DiffExp != "Non-significant"] <- processed_data$Gene[processed_data$DiffExp != "Non-significant"]
# 
# # Remove rows with NA in delabel
# labeled_data <- processed_data %>% na.omit(delabel)
# 
# # Sort data by log2ratio for up and down-regulated genes
# df_sorted_log2ratio_Down <- arrange(labeled_data, log2ratio)
# df_sorted_log2ratio_Up <- arrange(labeled_data, desc(log2ratio))
# 
# # Select top 30 down-regulated and top 50 up-regulated genes by log2ratio
# df_top30_log2ratio_Down <- slice(df_sorted_log2ratio_Down, 1:30)
# df_top50_log2ratio_Up <- slice(df_sorted_log2ratio_Up, 1:50)
# 
# # Sort data by FDR for up and down-regulated genes
# df_sorted_FDR_Down <- arrange(labeled_data, FDR)
# df_sorted_FDR_Up <- arrange(labeled_data, desc(FDR))
# 
# # Select top 50 genes by FDR for up and down-regulated genes
# df_top50_FDR_Down <- slice(df_sorted_FDR_Down, 1:50)
# df_top50_FDR_Up <- slice(df_sorted_FDR_Up, 1:50)

# Plot data with ggplot2
volcano_plot <- ggplot(processed_data, aes(log2ratio, -log10(FDR), color = DiffExp)) +
  geom_point() +
  scale_color_manual(values = c("blue", "red", "darkgrey"), name = NULL) +
  geom_vline(xintercept = c(-0.58, 0.58), col = "black", lwd = 0.5, lty = 2) +
  geom_hline(yintercept = -log10(0.05), col = "black", lwd = 0.5, lty = 2) +             
  # Add gene labels (uncomment the following lines if needed)
  # geom_text_repel(data = df_top30_log2ratio_Down, aes(label = delabel), color = "black", max.overlaps = 20) +
  # geom_text_repel(data = df_top50_log2ratio_Up, aes(label = delabel), color = "black", max.overlaps = 20) +
  # geom_text_repel(data = df_top50_FDR_Down, aes(label = delabel), color = "black", max.overlaps = 20) +
  # geom_text_repel(data = df_top50_FDR_Up, aes(label = delabel), color = "black", max.overlaps = 20) +
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

# Label genes based on differential expression
processed_data_5day$delabel <- NA
processed_data_5day$delabel[processed_data_5day$DiffExp != "Non-significant"] <- processed_data_5day$Gene[processed_data_5day$DiffExp != "Non-significant"]

# Remove rows with NA in delabel
labeled_data_5day <- processed_data_5day %>% na.omit(delabel)

# Sort data by log2ratio for up and down-regulated genes
df_sorted_log2ratio_Down_5day <- arrange(labeled_data_5day, log2ratio)
df_sorted_log2ratio_Up_5day <- arrange(labeled_data_5day, desc(log2ratio))

# Select top 50 down-regulated and up-regulated genes by log2ratio
df_top50_log2ratio_Down_5day <- slice(df_sorted_log2ratio_Down_5day, 1:50)
df_top50_log2ratio_Up_5day <- slice(df_sorted_log2ratio_Up_5day, 1:50)

# Sort data by FDR for down-regulated genes
df_sorted_FDR_Down_5day <- arrange(labeled_data_5day, FDR)

# Select top 50 genes by FDR for down-regulated genes
df_top50_FDR_Down_5day <- slice(df_sorted_FDR_Down_5day, 1:50)

# Plot data with ggplot2 for 5-day data
volcano_plot_5day <- ggplot(processed_data_5day, aes(log2ratio, -log10(FDR), color = DiffExp, label = delabel)) +
  geom_point() +
  scale_color_manual(values = c("blue", "red", "darkgrey"), name = NULL) +
  geom_vline(xintercept = c(-0.58, 0.58), col = "black", lwd = 0.5, lty = 2) +
  geom_hline(yintercept = -log10(0.05), col = "black", lwd = 0.5, lty = 2) +             
  # Add gene labels (uncomment the following lines if needed)
  # geom_text_repel(data = df_top50_log2ratio_Down_5day, aes(label = delabel), color = "black", max.overlaps = 20) +
  # geom_text_repel(data = df_top50_log2ratio_Up_5day, aes(label = delabel), color = "black", max.overlaps = 20) +
  # geom_text_repel(data = df_top50_FDR_Down_5day, aes(label = delabel), color = "black") +
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

