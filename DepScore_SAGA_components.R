# Set working directory
setwd("/Users/bgovinal/Documents/MDACC/Labmembers/Candy/DepMAP_ADA2B")

# Load necessary libraries
library(readxl)
library(tidyverse)
library(ggthemes)
library(ggpubr)
library(ggrepel)

# Read and preprocess data
df <- read_excel("SAGA dependency scores in MM.xlsx") %>% 
  rename(Cell_line = "Cell Line Name") %>% 
  gather(key = Gene_name, value = Expression, 2:22)

# Check the structure of the data
glimpse(df)

# Create a dataframe for gene module information
gene_module_info <- data.frame(
  Gene_name = c("TADA2B", "TADA3", "SGF29", "GCN5", "PCAF",
                "SUPT7L", "SUPT20H", "TADA1", "TAF5L", "TAF6L",   
                "ATXN7", "USP22", "ENY2", "ATXN7L3"),
  module = c("KAT", "KAT", "KAT", "KAT", "KAT",
             "Core", "Core", "Core", "Core", "Core",
             "DUB", "DUB", "DUB", "DUB")
)

# Merge data with gene module information
dat <- merge(df, gene_module_info, by = "Gene_name") %>%
  mutate(shape = case_when(Cell_line == "MM1S" ~ "MM1S",
                           Cell_line == "Other" ~ "Other MM cell lines",
                           TRUE ~ "Other MM cell lines"),
         Gene_name = recode(Gene_name, "TADA2B" = "ADA2B", "TADA3" = "ADA3", "TADA1" = "ADA1"))

dat$module <- factor(dat$module, levels = c("KAT", "Core", "DUB"))

# Check the head of the processed data
head(dat)

# Plotting the data
SAGA_effect <- dat %>% ggplot(aes(Gene_name, Expression)) +
  geom_boxplot(size = 0.3, outlier.colour = "white") +
  geom_jitter(aes(color = Cell_line), alpha = 0.7, position = position_jitter(0.1), show.legend = TRUE) +
  theme_par() +
  scale_size_manual(values = c(2, 1.5)) +
  scale_shape_manual(values = c(17, 16)) +
  scale_y_continuous(limits = c(-1.5, 0.1)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(y = "CRISPR Dependency Score",
       x = "",
       title = "",
       color = "MM Cell Lines") +
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 10),
        axis.ticks.length = unit(0.1, "cm"),
        axis.ticks = element_line(size = 0.2),
        axis.title.y = element_text(size = 10),
        plot.title = element_text(size = 10, face = "bold"),
        legend.position = "right",
        legend.text = element_text(size = 12)) +
  facet_wrap(~module, scales = "free")

# Enhance the plot
H_C_Dub <- SAGA_effect +
  guides(color = guide_legend(order = 1, override.aes = list(size = 4.2)),
         shape = guide_legend(order = 2),
         size = guide_legend(order = 3)) +
  theme(strip.text = element_text(face = "bold"))

# Display the plot
print(H_C_Dub)

# Save the plot
ggsave(H_C_Dub, filename = "SAGA_effect_box_Shape_Colored_Feb2024_ariel.tiff", width = 9, height = 5.5, dpi = 300)

### Supplemental Figures

# Read and preprocess supplemental data
df_sup <- read_excel("SAGA dependency scores in MM.xlsx") %>% 
  rename(Cell_line = "Cell Line Name") %>% 
  select(Cell_line, SF3B3:TAF12) %>% 
  gather(key = Gene_name, value = Expression, 2:8)

# Create a dataframe for supplemental gene module information
gene_module_info_sup <- data.frame(
  Gene_name = c("SF3B3", "SF3B5", "SUPT3H", "TAF9", "TAF10", "TAF12", "TRRAP"),
  module = c("Splicing", "Splicing", "Core", "Core", "Core", "Core", "TF-binding")
)

# Merge data with supplemental gene module information
dat_sup <- merge(df_sup, gene_module_info_sup, by = "Gene_name") %>%
  mutate(shape = case_when(Cell_line == "MM1S" ~ "MM1S",
                           Cell_line == "Other" ~ "Other MM cell lines",
                           TRUE ~ "Other MM cell lines"))

dat_sup$module <- factor(dat_sup$module, levels = c("Core", "Splicing", "TF-binding"))

# Check the head of the processed supplemental data
head(dat_sup)

# Plotting the supplemental data
SAGA_effect_sup <- dat_sup %>% ggplot(aes(Gene_name, Expression)) +
  geom_boxplot(size = 0.3, outlier.colour = "white") +
  geom_jitter(aes(color = Cell_line), alpha = 0.7, position = position_jitter(0.1), show.legend = TRUE) +
  theme_par() +
  scale_size_manual(values = c(2, 1.5)) +
  scale_shape_manual(values = c(17, 16)) +
  scale_y_continuous(limits = c(-3.5, 0.1)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(y = "CRISPR Dependency Score",
       x = "",
       title = "",
       color = "MM Cell Lines") +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 10),
        axis.ticks.length = unit(0.1, "cm"),
        axis.ticks = element_line(size = 0.2),
        axis.title.y = element_text(size = 10),
        plot.title = element_text(size = 10, face = "bold"),
        legend.position = "right",
        legend.text = element_text(size = 8)) +
  facet_wrap(~module, scales = "free")

# Enhance the supplemental plot
H_C_Dub_sup <- SAGA_effect_sup +
  guides(color = guide_legend(order = 1, override.aes = list(size = 4)),
         shape = guide_legend(order = 2),
         size = guide_legend(order = 3)) +
  theme(strip.text = element_text(face = "bold"),
        legend.direction = "horizontal",
        legend.position = "bottom")

# Display the supplemental plot
print(H_C_Dub_sup)

# Save the supplemental plot
ggsave(H_C_Dub_sup, filename = "SAGA_effect_box_Shape_Colored_Feb2024_supplement.tiff", width = 7, height = 6.5, dpi = 300)

### HAT Module

# Filter data for HAT components
HAT <- df %>% filter(Gene_name %in% c("TADA2B", "TADA3", "SGF29", "KAT2A"))

# Check the structure of the HAT data
glimpse(HAT)

# Plotting the HAT data
HAT_effect <- HAT %>% ggplot(aes(Gene_name, Expression)) +
  geom_boxplot(size = 0.2, outlier.colour = "white") +
  geom_jitter(aes(color = Cell_line), shape = 16, size = 1, position = position_jitter(0.1), show.legend = TRUE) +
  theme_par() +
  scale_y_continuous(limits = c(-1.5, 0.1)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(y = "Gene effect (Chronos) \nCRISPR (DepMap) Public 23Q2+Score",
       x = "",
       title = "HAT components") +
  theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 7),
        axis.ticks.length = unit(0.1, "cm"),
        axis.ticks = element_line(size = 0.2),
        axis.title.y = element_text(size = 7, face = "bold"),
        plot.title = element_text(size = 10))

# Display the HAT plot
print(HAT_effect)

# Save the HAT plot
ggsave(HAT_effect, filename = "HAT_effect.tiff", width = 3, height = 3, dpi = 300)
