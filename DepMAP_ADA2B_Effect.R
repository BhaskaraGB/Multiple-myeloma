# Clear the workspace
rm(list = ls())

# Set working directory
setwd("/Users/bgovinal/Documents/MDACC/Labmembers/Candy/DepMAP_ADA2B")

# Load necessary libraries
library(readxl)
library(tidyverse)
library(ggthemes)
library(ggpubr)
library(ggrepel)

# Read and preprocess data
df1 <- read.csv("DepMAP_ADA2B.csv") %>% 
  rename(Expression = "Expression.Public.23Q2") %>% 
  select(Cell_Line, Lineage, CRISPR_Score, Expression) %>% 
  na.omit() %>% 
  group_by(Lineage) %>% 
  summarise(mean_CrisprScore = mean(CRISPR_Score),
            mean_Expression = mean(Expression)) %>% 
  filter(mean_Expression > 1) %>% 
  mutate(color = case_when(Lineage == "Plasma_cells" ~ "blue",
                           Lineage == "Other" ~ "gray85",
                           TRUE ~ "black"),
         status = ifelse(mean_CrisprScore <= -0.42, "low", "high"))

# Glimpse of the processed data
glimpse(df1)

# Label differentially expressed genes
df1 <- df1 %>%
  mutate(label = ifelse(status != "high", Lineage, NA))

# Plotting the data
TADA2B_effect <- df1 %>%
  ggplot(aes(mean_CrisprScore, mean_Expression, fill = Lineage)) +
  geom_point(alpha = 0.7, shape = 21, size = 5, show.legend = TRUE) +
  scale_y_continuous(limits = c(2, 4), breaks = seq(2, 4, 1)) +
  scale_x_continuous(limits = c(-1, 0.5)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  geom_text_repel(aes(label = label), size = 4, hjust = 1.5, vjust = 0) +
  labs(y = "ADA2B Expression \n(Log2(TPM+1))",
       x = "ADA2B Effect \nCRISPR Dependency Score") +
  theme_par() +
  theme(axis.title = element_text(face = "bold"),
        axis.text = element_text(),
        axis.ticks.length = unit(0.1, "cm"),
        legend.title = element_text(face = "bold"))

# Display the plot
print(TADA2B_effect)

# Save the plot
ggsave(TADA2B_effect, filename = "TADA2B_effect_with_leukemia.tiff", width = 7, height = 7, dpi = 300)
