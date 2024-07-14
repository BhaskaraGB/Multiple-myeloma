Here is a README file for your GitHub repository:
  
  ### README.md
  
  ```markdown
# ATAC-Seq Analysis of Gene Expression Changes

This repository contains R scripts and data for analyzing ATAC-seq peaks associated with upregulated and downregulated genes in a specific experiment. The analysis involves reading peak information, normalizing read counts to counts per million (CPM), and generating plots to visualize the data.

## Files

- `peakinfo.csv`: Contains peak information with columns for chromosome, start, end, and associated genes.
- `number_reads.csv`: Contains read counts for various samples.
- `NT-3d-Dox_vs_s2B-3d-Dox_edgeR_FDR0.05.final.csv`: Contains differential expression analysis results with log2 ratios for 3-day treatment.
- `NT-5d-Dox_vs_s2B-5d-Dox_edgeR_FDR0.05.final.csv`: Contains differential expression analysis results with log2 ratios for 5-day treatment.
- `2B-3d-Dox-1.over_peak_center.txt`, `2B-3d-Dox-2.over_peak_center.txt`, `NT-3d-Dox-1.over_peak_center.txt`, `NT-3d-Dox-2.over_peak_center.txt`: Read coverage data for various samples (3-day treatment).
- `2B-5d-Dox-1.over_peak_center.txt`, `2B-5d-Dox-2.over_peak_center.txt`, `NT-5d-Dox-1.over_peak_center.txt`, `NT-5d-Dox-2.over_peak_center.txt`: Read coverage data for various samples (5-day treatment).

## Script

The R script reads data from the provided files, processes the data to identify and normalize upregulated and downregulated genes, and generates plots to visualize the normalized read coverage.

### Steps

1. **Set Working Directory and Load Libraries**:
  ```r
rm(list = ls())
setwd("/path/to/your/directory")
library(tidyverse)
library(ggplot2)
library(readxl)
library(ggthemes)
library(patchwork)
```

2. **Read Peak Information**:
  ```r
PeakList <- read.csv("peakinfo.csv") %>% 
  mutate(PeakID = paste(Chr, Start, End, sep = "_"))
```

3. **Read Coverage Information**:
  ```r
Coverage <- read.csv("number_reads.csv")
```

4. **Read and Categorize Differential Expression Data (3-Day Treatment)**:
  ```r
Exp3D_NTvs2B_Dox <- read.csv("NT_3D_D0x_vs_s2B_3d_Dox_edgeR_FDR0.05_Final_categorynames.csv") %>% 
  mutate(category = ifelse(log2ratio >= 0.585, "Up", ifelse(log2ratio <= -0.585, "Down", NA)))
table(Exp3D_NTvs2B_Dox$category)
```

5. **Filter Peaks Based on Differential Expression (3-Day Treatment)**:
  ```r
PeakListExp3D <- PeakList[which(PeakList$Gene.3d %in% Exp3D_NTvs2B_Dox$gene_name),]
```

6. **Read and Process Samples for Downregulated Expression (3-Day Treatment)**:
  ```r
D3_2B_REP1 <- read.delim("2B-3d-Dox-1.over_peak_center.txt", header = F, sep = "\t") %>% 
  mutate(PeakID = PeakList$PeakID) %>% 
  filter(PeakID %in% PeakList$PeakID[-which(is.na(PeakList$Gene.3d))]) %>% 
  filter(PeakID %in% PeakListExp3D$PeakID)

PeakOrder <- as.data.frame(cbind(PeakID = D3_2B_REP1$PeakID)) %>% 
  left_join(., PeakListExp3D[, c('Gene.3d', 'PeakID')]) %>% 
  left_join(., Exp3D_NTvs2B_Dox[, c('gene_name', 'category')], by = c("Gene.3d" = "gene_name"))

D3_2B_REP2 <- read.delim("2B-3d-Dox-2.over_peak_center.txt", header = F, sep = "\t") %>% 
  mutate(PeakID = PeakList$PeakID) %>% 
  filter(PeakID %in% PeakList$PeakID[-which(is.na(PeakList$Gene.3d))]) %>% 
  filter(PeakID %in% PeakListExp3D$PeakID) %>% 
  arrange(factor(PeakID, levels = PeakOrder$PeakID))
```

7. **Normalize Read Coverage for Downregulated Expression (3-Day Treatment)**:
  ```r
dat3Day_2B <- as.data.frame(cbind(
  D3_2B_REP1 = colMeans(D3_2B_REP1[which(PeakOrder$category == "Down"), 1:500]) * 1000000 / (Coverage$TotalReads[which(Coverage$Sample == "2B-3d-Dox-1")]),
  D3_2B_REP2 = colMeans(D3_2B_REP2[which(PeakOrder$category == "Down"), 1:500]) * 1000000 / (Coverage$TotalReads[which(Coverage$Sample == "2B-3d-Dox-2")])
)) %>% 
  as_tibble() %>% 
  mutate(MEAN = rowMeans(., na.rm = T)) %>% 
  mutate(Bin = seq(1:500)) %>% 
  mutate(Sample = "shADA2B") %>% 
  mutate(Time = "Day-3") %>% 
  select(MEAN, Bin, Sample, Time)
```

8. **Generate Plot for Downregulated Expression (3-Day Treatment)**:
  ```r
DownExp_ATAC_Day3 <- as_tibble(rbind(dat3Day_2B, dat3Day_NT)) %>% 
  ggplot(aes(x = Bin, y = MEAN)) +
  scale_x_continuous(breaks = c(0, 500), labels = c("-5kb", "+5kb")) +
  geom_line(aes(color = Sample), size = 1.4) +
  labs(x = "", y = "Normalized Read Coverage", title = "ATAC-seq Peak of \nDownregulated Genes\n (3d, 533)") +
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

DownExp_ATAC_Day3
```

9. **Repeat Steps 6-8 for Upregulated Expression (3-Day Treatment)**
  
  10. **Repeat Steps 4-9 for 5-Day Treatment**
  
  11. **Combine and Save Plots**:
  ```r
Exp_ATAC_day3 <- DownExp_ATAC_Day3 + UpExp_ATAC_Day3
ggsave(Exp_ATAC_day3, filename = "ATAC_Exp_log2ratio_day3_July14_24.tiff", width = 12, height = 7, dpi = 300)

Exp_ATAC_D5 <- DownExp_ATAC_D5 + UpExp_ATAC_D5
ggsave(Exp_ATAC_D5, filename = "ATAC_Exp_log2ratio_DAY5_July14_24.tiff", width = 12, height = 7, dpi = 300)
```

## Dependencies

- `tidyverse`
- `ggplot2`
- `readxl`
- `ggthemes`
- `patchwork`

## How to Run the Script

1. Ensure all the required input files are in the same directory as the script.
2. Load the necessary R libraries.
3. Execute the script step by step in an R environment or run it as an R script.
4. The script will generate plots visualizing the normalized read coverage for upregulated and downregulated genes.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

```

This README file provides a comprehensive and nicely formatted overview of your analysis script, including a description of the files, steps to run the