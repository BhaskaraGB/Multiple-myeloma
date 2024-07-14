# Clear the workspace
rm(list = ls())

# Set the working directory
setwd("/Users/bgovinal/Documents/MDACC/Labmembers/Candy/Metaplots_ATAC/OneDrive_1_12-6-2023/08302022_ATAC ADA2B MM")

# Load necessary libraries
library(tidyverse)
library(ggplot2)
library(readxl)
library(ggthemes)
library(patchwork)

# Read peak information and create a unique PeakID
PeakList <- read.csv("peakinfo.csv") %>% 
  mutate(PeakID = paste(Chr, Start, End, sep = "_"))

# Read coverage information
Coverage <- read.csv("number_reads.csv")

### 5-Day Treatment

# Read and categorize differential expression data (5-Day)
Exp5D_NTvs2B_Dox <- read.csv("NT-5d-Dox_vs_s2B-5d-Dox_edgeR_FDR0.05.final.csv") %>% 
  mutate(category = ifelse(log2ratio >= 0.585, "Up", ifelse(log2ratio <= -0.585, "Down", NA)))

table(Exp5D_NTvs2B_Dox$category)

# Filter peaks based on differential expression data
PeakListExp5D <- PeakList[which(PeakList$Gene.3d %in% Exp5D_NTvs2B_Dox$gene_name),]

# Downregulated expression - Day 5
# 2B
D5_2B_REP1 <- read.delim("2B-5d-Dox-1.over_peak_center.txt", header = F, sep = "\t") %>% 
  mutate(PeakID = PeakList$PeakID) %>% 
  filter(PeakID %in% PeakList$PeakID[-which(is.na(PeakList$Gene.5d))]) %>% 
  filter(PeakID %in% PeakListExp5D$PeakID)

PeakOrder <- as.data.frame(cbind(PeakID = D5_2B_REP1$PeakID)) %>% 
  left_join(., PeakListExp5D[, c('Gene.5d', 'PeakID')]) %>% 
  left_join(., Exp5D_NTvs2B_Dox[, c('gene_name', 'category')], by = c("Gene.5d" = "gene_name"))

D5_2B_REP2 <- read.delim("2B-5d-Dox-2.over_peak_center.txt", header = F, sep = "\t") %>% 
  mutate(PeakID = PeakList$PeakID) %>% 
  filter(PeakID %in% PeakList$PeakID[-which(is.na(PeakList$Gene.5d))]) %>% 
  filter(PeakID %in% PeakListExp5D$PeakID) %>% 
  arrange(factor(PeakID, levels = PeakOrder$PeakID))

# Select "Down" category and normalize read coverage
dat5Day_2B <- as.data.frame(cbind(
  D5_2B_REP1 = colMeans(D5_2B_REP1[which(PeakOrder$category == "Down"), 1:500]) * 1000000 / (Coverage$TotalReads[which(Coverage$Sample == "2B-5d-Dox-1")]),
  D5_2B_REP2 = colMeans(D5_2B_REP2[which(PeakOrder$category == "Down"), 1:500]) * 1000000 / (Coverage$TotalReads[which(Coverage$Sample == "2B-5d-Dox-2")])
)) %>% 
  as_tibble() %>% 
  mutate(MEAN = rowMeans(., na.rm = T)) %>% 
  mutate(Bin = seq(1:500)) %>% 
  mutate(Sample = "shADA2B") %>% 
  mutate(Time = "Day-5") %>% 
  select(MEAN, Bin, Sample, Time)

# NT
D5_NT_REP1 <- read.delim("NT-5d-Dox-1.over_peak_center.txt", header = F, sep = "\t") %>% 
  mutate(PeakID = PeakList$PeakID) %>% 
  filter(PeakID %in% PeakList$PeakID[-which(is.na(PeakList$Gene.5d))]) %>% 
  filter(PeakID %in% PeakListExp5D$PeakID) %>% 
  arrange(factor(PeakID, levels = PeakOrder$PeakID))

D5_NT_REP2 <- read.delim("NT-5d-Dox-2.over_peak_center.txt", header = F, sep = "\t") %>% 
  mutate(PeakID = PeakList$PeakID) %>% 
  filter(PeakID %in% PeakList$PeakID[-which(is.na(PeakList$Gene.5d))]) %>% 
  filter(PeakID %in% PeakListExp5D$PeakID) %>% 
  arrange(factor(PeakID, levels = PeakOrder$PeakID))

dat5Day_NT <- as.data.frame(cbind(
  D5_NT_REP1 = colMeans(D5_NT_REP1[which(PeakOrder$category == "Down"), 1:500]) * 1000000 / (Coverage$TotalReads[which(Coverage$Sample == "NT-5d-Dox-1")]),
  D5_NT_REP2 = colMeans(D5_NT_REP2[which(PeakOrder$category == "Down"), 1:500]) * 1000000 / (Coverage$TotalReads[which(Coverage$Sample == "NT-5d-Dox-2")])
)) %>% 
  as_tibble() %>% 
  mutate(MEAN = rowMeans(., na.rm = T)) %>% 
  mutate(Bin = seq(1:500)) %>% 
  mutate(Sample = "shNT") %>% 
  mutate(Time = "Day-5") %>% 
  select(MEAN, Bin, Sample, Time)

# Plot for Downregulated Expression (5-Day)
DownExp_ATAC_D5 <- as_tibble(rbind(dat5Day_2B, dat5Day_NT)) %>% 
  ggplot(aes(x = Bin, y = MEAN)) +
  scale_x_continuous(breaks = c(0, 500), labels = c("-5kb", "+5kb")) +
  geom_line(aes(color = Sample), size = 1.4) +
  labs(x = "", y = "Normalized Read Coverage", title = "ATAC-seq peak of \n Downregulated Genes \n (5d, 474)") +
  scale_color_manual(values = c("blue", "red")) +
  scale_y_continuous(limits = c(0, 2), breaks = seq(0, 2, 0.5)) +
  theme_par() +
  theme(
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 18),
    axis.title.x = element_blank(),
    axis.text.x = element_text(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 16, vjust = 1, color = "black"),
    legend.position = c(0.75, 0.88),
    legend.title = element_blank(),
    legend.text = element_text(size = 12, vjust = 1, color = "black", face = "bold")
  )

DownExp_ATAC_D5

#----------------------------------------------------------------------------------

### Upregulated expression - Day 5

# 2B
D5_2B_REP1 <- read.delim("2B-5d-Dox-1.over_peak_center.txt", header = F, sep = "\t") %>% 
  mutate(PeakID = PeakList$PeakID) %>% 
  filter(PeakID %in% PeakList$PeakID[-which(is.na(PeakList$Gene.5d))]) %>% 
  filter(PeakID %in% PeakListExp5D$PeakID)

PeakOrder <- as.data.frame(cbind(PeakID = D5_2B_REP1$PeakID)) %>% 
  left_join(., PeakListExp5D[, c('Gene.5d', 'PeakID')]) %>% 
  left_join(., Exp5D_NTvs2B_Dox[, c('gene_name', 'category')], by = c("Gene.5d" = "gene_name"))

D5_2B_REP2 <- read.delim("2B-5d-Dox-2.over_peak_center.txt", header = F, sep = "\t") %>% 
  mutate(PeakID = PeakList$PeakID) %>% 
  filter(PeakID %in% PeakList$PeakID[-which(is.na(PeakList$Gene.5d))]) %>% 
  filter(PeakID %in% PeakListExp5D$PeakID) %>% 
  arrange(factor(PeakID, levels = PeakOrder$PeakID))

# Select "Up" category and normalize read coverage
dat5Day_2B <- as.data.frame(cbind(
  D5_2B_REP1 = colMeans(D5_2B_REP1[which(PeakOrder$category == "Up"), 1:500]) * 1000000 / (Coverage$TotalReads[which(Coverage$Sample == "2B-5d-Dox-1")]),
  D5_2B_REP2 = colMeans(D5_2B_REP2[which(PeakOrder$category == "Up"), 1:500]) * 1000000 / (Coverage$TotalReads[which(Coverage$Sample == "2B-5d-Dox-2")])
)) %>% 
  as_tibble() %>% 
  mutate(MEAN = rowMeans(., na.rm = T)) %>% 
  mutate(Bin = seq(1:500)) %>% 
  mutate(Sample = "shADA2B") %>% 
  mutate(Time = "Day-5") %>% 
  select(MEAN, Bin, Sample, Time)

# NT
D5_NT_REP1 <- read.delim("NT-5d-Dox-1.over_peak_center.txt", header = F, sep = "\t") %>% 
  mutate(PeakID = PeakList$PeakID) %>% 
  filter(PeakID %in% PeakList$PeakID[-which(is.na(PeakList$Gene.5d))]) %>% 
  filter(PeakID %in% PeakListExp5D$PeakID) %>% 
  arrange(factor(PeakID, levels = PeakOrder$PeakID))

D5_NT_REP2 <- read.delim("NT-5d-Dox-2.over_peak_center.txt", header = F, sep = "\t") %>% 
  mutate(PeakID = PeakList$PeakID) %>% 
  filter(PeakID %in% PeakList$PeakID[-which(is.na(PeakList$Gene.5d))]) %>% 
  filter(PeakID %in% PeakListExp5D$PeakID) %>% 
  arrange(factor(PeakID, levels = PeakOrder$PeakID))

dat5Day_NT <- as.data.frame(cbind(
  D5_NT_REP1 = colMeans(D5_NT_REP1[which(PeakOrder$category == "Up"), 1:500]) * 1000000 / (Coverage$TotalReads[which(Coverage$Sample == "NT-5d-Dox-1")]),
  D5_NT_REP2 = colMeans(D5_NT_REP2[which(PeakOrder$category == "Up"), 1:500]) * 1000000 / (Coverage$TotalReads[which(Coverage$Sample == "NT-5d-Dox-2")])
)) %>% 
  as_tibble() %>% 
  mutate(MEAN = rowMeans(., na.rm = T)) %>% 
  mutate(Bin = seq(1:500)) %>% 
  mutate(Sample = "shNT") %>% 
  mutate(Time = "Day-5") %>% 
  select(MEAN, Bin, Sample, Time)

# Plot for Upregulated Expression (5-Day)
UpExp_ATAC_D5 <- as_tibble(rbind(dat5Day_2B, dat5Day_NT)) %>% 
  ggplot(aes(x = Bin, y = MEAN)) +
  scale_x_continuous(breaks = c(0, 500), labels = c("-5kb", "+5kb")) +
  geom_line(aes(color = Sample), size = 1.4) +
  labs(x = "", y = "Normalized Read Coverage", title = "ATAC-seq peak of \n Upregulated Genes \n (5d, 502)") +
  scale_color_manual(values = c("blue", "red")) +
  theme_par() +
  theme(
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 18),
    axis.title.x = element_blank(),
    axis.text.x = element_text(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 16, vjust = 1, color = "black"),
    legend.position = c(0.75, 0.88),
    legend.title = element_blank(),
    legend.text = element_text(size = 12, vjust = 1, color = "black", face = "bold")
  )

UpExp_ATAC_D5

Exp_ATAC_D5 <- DownExp_ATAC_D5 + UpExp_ATAC_D5
ggsave(Exp_ATAC_D5, filename = "ATAC_Exp_log2ratio_DAY5_July14_24.tiff", width = 12, height = 7, dpi = 300)


