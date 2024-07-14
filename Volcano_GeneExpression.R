# Clear the workspace
rm(list = ls())

# Set working directory
setwd("/Users/bgovinal/Documents/MDACC/Labmembers/Candy/gene_expression")

# Load necessary libraries
library(ggplot2)
library(ggthemes)
library(ggrepel)
library(tidyverse)

### 3-Day Treatment Analysis

# Read input CSV file for 3-day data
data_3day <- read.csv("NT_3D_D0x_vs_s2B_3d_Dox_edgeR_All_Final_categorynames.csv")

# Process data: Select relevant columns, categorize differential expression, and label target genes
processed_data_3day <- data_3day %>% 
  select(gene_name, log2ratio, FDR) %>% 
  mutate(DiffExp = ifelse(FDR <= 0.05 & abs(log2ratio) >= 0, 
                          ifelse(log2ratio >= 0, "Increased", "Decreased"),
                          "Non-significant")) %>% 
  mutate(gene_name = recode(gene_name, "TADA2B" = "ADA2B")) %>% 
  mutate(label = ifelse(gene_name == "ADA2B", "ADA2B", NA)) %>% 
  mutate(DiffExp = ifelse(gene_name == "ADA2B", "Target", DiffExp))

# Display the number of differentially expressed genes
table(processed_data_3day$DiffExp)

# Plot data for 3-day treatment
plot_3day <- processed_data_3day %>% ggplot(aes(log2ratio, -log10(FDR), color = DiffExp, label = label)) +
  geom_point() +
  geom_vline(xintercept = c(-0.58, 0.58), col = "black", lwd = 0.5, lty = 2) +
  geom_hline(yintercept = -log10(0.05), col = "black", lwd = 0.5, lty = 2) +
  scale_color_manual(values = c("Decreased" = "blue", "Increased" = "red", 
                                "Non-significant" = "darkgrey", "Target" = "skyblue"),
                     name = NULL) +
  scale_x_continuous(breaks = seq(-6, 6, 2)) +
  scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, 20)) +
  labs(title = "NT-3Day-Dox vs ADA2b",
       x = expression("Log"[2] ~ "FC"),
       y = expression("Log"[10] ~ "FDR")) +
  theme_par() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        title = element_text(size = 20),
        legend.text = element_text(size = 15)) +
  guides(size = FALSE, color = guide_legend(override.aes = list(size = 5)))

# Display the plot for 3-day treatment
print(plot_3day)

### 5-Day Treatment Analysis

# Read input CSV file for 5-day data
data_5day <- read.csv("NT-5d-Dox_vs_s2B-5d-Dox_edgeR_All.final.csv")

# Process data: Select relevant columns, categorize differential expression, and label target genes
processed_data_5day <- data_5day %>% 
  select(gene_name, log2ratio, FDR) %>% 
  mutate(DiffExp = ifelse(FDR <= 0.05 & abs(log2ratio) >= 0, 
                          ifelse(log2ratio >= 0, "Increased", "Decreased"),
                          "Non-significant")) %>% 
  mutate(gene_name = recode(gene_name, "TADA2B" = "ADA2B")) %>% 
  mutate(label = ifelse(gene_name == "ADA2B", "ADA2B", NA)) %>% 
  mutate(DiffExp = ifelse(gene_name == "ADA2B", "Target", DiffExp))

# Plot data for 5-day treatment
plot_5day <- processed_data_5day %>% ggplot(aes(log2ratio, -log10(FDR), color = DiffExp, label = label)) +
  geom_point() +
  geom_vline(xintercept = c(-0.58, 0.58), col = "black", lwd = 0.5, lty = 2) +
  geom_hline(yintercept = -log10(0.05), col = "black", lwd = 0.5, lty = 2) +
  scale_color_manual(values = c("Decreased" = "blue", "Increased" = "red", 
                                "Non-significant" = "darkgrey", "Target" = "skyblue"),
                     name = NULL) +
  scale_x_continuous(breaks = seq(-6, 6, 2)) +
  scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, 20)) +
  labs(title = "2B vs NT 5d-Dox",
       x = expression("Log"[2] ~ "FC"),
       y = expression("Log"[10] ~ "FDR")) +
  theme_par() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        title = element_text(size = 20),
        legend.text = element_text(size = 15)) +
  guides(size = FALSE, color = guide_legend(override.aes = list(size = 5)))

# Display the plot for 5-day treatment
print(plot_5day)

# Combine and display the plots for 3-day and 5-day treatments
combined_volcano <- plot_3day + plot_5day
print(combined_volcano)

# Save the combined plot
ggsave(combined_volcano, file = "Volcano_expression_v3_nonames_feb13_2024.tiff", dpi = 300, width = 18, height = 8)
