# Set working directory
setwd("/Users/bgovinal/Documents/MDACC/Labmembers/Candy/MAF_targets_heatmap/GB_analysis/")

# Load necessary libraries
library(pheatmap)
library(tidyverse)

# Load MAF target data and expression data for 3-day and 5-day Dox treatments
df1 <- read.csv("MAF_targets.csv")
df2 <- read.csv("NT_s2B_3d_Dox_edgeR_All.final.csv")
df3 <- read.csv("NT_s2B_5d_Dox_edgeR_All.final.csv")

# Merge MAF target data with 3-day Dox expression data and save to CSV
merge_df1_df2 <- merge(df1, df2, by = "gene_name")
write.csv(merge_df1_df2, file = "MAF_targets_exp_3d_Dox.csv")

# Merge MAF target data with 5-day Dox expression data and save to CSV
merge_df1_df3 <- merge(df1, df3, by = "gene_name")
MAF_target_exp_unique <- merge_df1_df3 %>% distinct(gene_name, .keep_all = TRUE)
write.csv(MAF_target_exp_unique, file = "MAF_targets_exp_5d_Dox.csv")

# Load merged data for 3-day and 5-day Dox treatments
df4 <- read.csv("MAF_targets_exp_3d_Dox.csv") %>% 
  select(gene_name, exp1, exp2, log2ratio) %>% 
  rename(exp1_3d = exp1, exp2_3d = exp2, log2ratio_3d = log2ratio)

df5 <- read.csv("MAF_targets_exp_5d_Dox.csv") %>% 
  select(gene_name, exp1, exp2, log2ratio) %>% 
  rename(exp1_5d = exp1, exp2_5d = exp2, log2ratio_5d = log2ratio)

# Merge 3-day and 5-day expression data
MAF_target_3d_5d_Dox <- merge(df4, df5, by = "gene_name")
write.csv(MAF_target_3d_5d_Dox, file = "MAF_target_3d_5d_Dox.csv", row.names = FALSE)

# Load ADA2B binding data and merge with the combined expression data
df6 <- read.csv("Flag_3xHA-ADA2B_HA_TechRep-2_q0.05.geneannotation.csv") %>% 
  select(gene_name) %>% 
  distinct()

df7 <- read.csv("MAF_target_3d_5d_Dox.csv")
ada2bind <- merge(df6, df7, by = "gene_name")
write.csv(ada2bind, file = "ada2b_bind_gene.csv")

# Prepare data for heatmap
mat <- read.csv("MAF_target_3d_5d_Dox.csv")
mat_exp <- mat %>% 
  select(gene_name, log2ratio_3d, log2ratio_5d, type) %>% 
  rename(day3 = log2ratio_3d, day5 = log2ratio_5d)

# Generate heatmap
heatPlot <- mat_exp %>% 
  gather(key = day, value = Expression, day3:day5) %>% 
  ggplot(aes(y = gene_name, x = day)) +
  geom_tile(aes(fill = Expression)) +
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  geom_segment(aes(x = 2.53, xend = 2.53, y = 0.5, yend = 20.5), color = "seagreen", linewidth = 3) +
  geom_segment(aes(x = 2.53, xend = 2.53, y = 20.5, yend = 36.5), color = "black", linewidth = 3) +
  labs(y = "", x = "Days of Dox Treatment", fill = "Log2 ratio") +
  theme_void() +
  theme(axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8, face = "bold",
                                   colour = c(rep("seagreen", 20), rep("black", 16))))

# Save the heatmap
ggsave(heatPlot, filename = "heatplot_v5_Feb12_2024.pdf", width = 7, height = 6, dpi = 300)
