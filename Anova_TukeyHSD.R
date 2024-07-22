rm(list = ls())

# Load required libraries
library(tidyverse)
library(readxl)
library(ggthemes)
library(broom)
library(dplyr)
library(multcompView)

# Set working directory
setwd("/Users/bgovinal/Documents/MDACC/Labmembers/Candy/Resuce_expt/July2024/")


perform_anova_tukey <- function(data, day, treat_condition) {
  # Filter the data based on cell type and treatment condition
  filtered_data <- data %>%
    filter(Day == day & Treat == treat_condition)
  
  # Perform ANOVA
  res.aov <- aov(Cell_number ~ Genotype, data = filtered_data)
  aov_summary <- tidy(res.aov)
  
  # Perform Tukey HSD test
  tukey_results <- TukeyHSD(res.aov)
  
  list(summary = aov_summary, tukey = tukey_results)
}


# Load and glimpse MM1S cell viability data
MM1S_viability <- read_excel("Fig 1E_Cell_viability.xlsx", sheet = "1st experiment") %>% 
  select(Genotype, Treat, Cell_number, Day) %>% 
  drop_na()

glimpse(MM1S_viability)

# Function to perform ANOVA and Tukey HSD on 9day, No Dox
anova_day9_No_Dox <- perform_anova_tukey(MM1S_viability, "day9", "No_Dox") 
print(anova_day9_No_Dox$summary)
print(anova_day9_No_Dox$tukey)

# Function to perform ANOVA and Tukey HSD on 9day, +Dox
anova_day9_Dox <- perform_anova_tukey(MM1S_viability, "day9", "Dox") 
print(anova_day9_Dox$summary)
print(anova_day9_Dox$tukey)

# Function to perform ANOVA and Tukey HSD on 12day, No Dox
anova_day12_No_Dox <- perform_anova_tukey(MM1S_viability, "day12", "No_Dox") 
print(anova_day12_No_Dox$summary)
print(anova_day12_No_Dox$tukey)

# Function to perform ANOVA and Tukey HSD on 9day, +Dox
anova_day12_Dox <- perform_anova_tukey(MM1S_viability, "day12", "Dox") 
print(anova_day12_Dox$summary)
print(anova_day12_Dox$tukey)



# Load and glimpse cell cycle data
Cell_cycle <- read_excel("Fig 3G cell cycle.xlsx") %>% 
  select(Genotype, Cell_cycle, Value)

glimpse(Cell_cycle)

# Function to perform ANOVA and Tukey HSD for cell cycle phases
perform_anova <- function(data, phase) {
  phase_data <- data %>% filter(Cell_cycle == phase)
  phase_aov <- aov(Value ~ Genotype, data = phase_data)
  phase_aov_summary <- tidy(phase_aov)
  phase_tukey <- TukeyHSD(phase_aov)
  
  list(summary = phase_aov_summary, tukey = phase_tukey)
}

# ANOVA and Tukey HSD for G0/G1 phase
G0_G1_results <- perform_anova(Cell_cycle, "G0/G1")
print(G0_G1_results$summary)
print(G0_G1_results$tukey)

# ANOVA and Tukey HSD for S phase
S_phase_results <- perform_anova(Cell_cycle, "S")
print(S_phase_results$summary)
print(S_phase_results$tukey)

# ANOVA and Tukey HSD for G2/M phase
G2_M_phase_results <- perform_anova(Cell_cycle, "G2/M")
print(G2_M_phase_results$summary)
print(G2_M_phase_results$tukey)


#Fig4c, S6C CPTH6
#------------------------------------------------------------------------------

CPTH6_10uM <- read_excel("Fig 4C, S6C CPTH6.xlsx", sheet = "10uM") %>% 
  select(1:4) %>% 
  drop_na()

glimpse(CPTH6_10uM)

perform_anova_tukey <- function(data, day) {
  day_data <- data %>% filter(Day == day)
  res_aov <- aov(Values ~ Treat, data = day_data)
  res_aov_summary <- tidy(res_aov)
  tukey_results <- TukeyHSD(res_aov)
  
  list(summary = res_aov_summary, tukey = tukey_results)
}

# Perform ANOVA and Tukey HSD for 3-day data
DAY3_results <- perform_anova_tukey(CPTH6_10uM, "Day3")

# Print the results for 3-day data
print(DAY3_results$summary)
print(DAY3_results$tukey)

# Perform ANOVA and Tukey HSD for 6-day data
DAY6_results <- perform_anova_tukey(CPTH6_10uM, "Day6")

# Print the results for 6-day data
print(DAY6_results$summary)
print(DAY6_results$tukey)


#####Fig S6C CPTH6
###----------------------------------------------------------------------------
#Fig4c, S6C CPTH6
#------------------------------------------------------------------------------

rm(list = ls())

##load data
CPTH6_1_2_5uM <- read_excel("Fig 4C, S6C CPTH6.xlsx", sheet = "1,2,5uM") %>% 
  select(Treat, Day, Values) %>% 
  drop_na()

glimpse(CPTH6_1_2_5uM)

##Create function to perform Anova and TukeyHSD

perform_anova_tukey <- function(data, day) {
  day_data <- data %>% filter(Day == day)
  res_aov <- aov(Values ~ Treat, data = day_data)
  res_aov_summary <- tidy(res_aov)
  tukey_results <- TukeyHSD(res_aov)
  
  list(summary = res_aov_summary, tukey = tukey_results)
}


# Perform ANOVA and Tukey HSD for 3-day data
DAY3_results <- perform_anova_tukey(CPTH6_1_2_5uM,  "Day3") 

# Print the results for 3-day data
print(DAY3_results$summary)
print(DAY3_results$tukey)

# Perform ANOVA and Tukey HSD for 6-day data
DAY6_results <- perform_anova_tukey(CPTH6_1_2_5uM,  "Day6")

# Print the results for 6-day data
print(DAY6_results$summary)
print(DAY6_results$tukey)


# Fig 6G FHA rescue
#-------------------------------------------------------------------------------
FHA_rescue <- read_excel("Fig 6G FHA rescue.xlsx", sheet = "Rescue") 

# Prepare data
df_rescue <- FHA_rescue %>% 
  select(Genotypes, Cell_count) %>% 
  drop_na()

df_rescue$Cell_count <- as.numeric(df_rescue$Cell_count)

glimpse(df_rescue)

# ANOVA for Rescue Experiment
rescue_aov <- aov(Cell_count ~ Genotypes, data = df_rescue)
rescue_aov_summary <- tidy(rescue_aov)
rescue_tukey <- TukeyHSD(rescue_aov)

# Print ANOVA summary and Tukey HSD results
print(rescue_aov_summary)
print(rescue_tukey)


## Fig S1D Crispr KO
#-------------------------------------------------------------------------------

CRISPR <-  read_excel("Fig S1D Crispr KO.xlsx") %>% 
  select(1:4) %>% 
  drop_na()

# Remove the space in the Days column
CRISPR$Day <- gsub("Day\\s+(\\d+)", "Day\\1", CRISPR$Day)

glimpse(CRISPR)

perform_anova_tukey <- function(data, day) {
  day_data <- data %>% filter(Day == day)
  res_aov <- aov(Values ~ Genotype, data = day_data)
  res_aov_summary <- tidy(res_aov)
  tukey_results <- TukeyHSD(res_aov)
  
  list(summary = res_aov_summary, tukey = tukey_results)
}


# Perform ANOVA and Tukey HSD for 0-day data
DAY0_results <- perform_anova_tukey(CRISPR, "Day0")

# Print the results for 0-day data
print(DAY0_results$summary)
print(DAY0_results$tukey)


# Perform ANOVA and Tukey HSD for 2-day data
DAY2_results <- perform_anova_tukey(CRISPR, "Day2")

# Print the results for 2-day data
print(DAY2_results$summary)
print(DAY2_results$tukey)

# Perform ANOVA and Tukey HSD for 2-day data
DAY5_results <- perform_anova_tukey(CRISPR, "Day5")

# Print the results for 2-day data
print(DAY5_results$summary)
print(DAY5_results$tukey)


# Perform ANOVA and Tukey HSD for 2-day data
DAY8_results <- perform_anova_tukey(CRISPR, "Day8")

# Print the results for 2-day data
print(DAY8_results$summary)
print(DAY8_results$tukey)


# Perform ANOVA and Tukey HSD for 2-day data
DAY10_results <- perform_anova_tukey(CRISPR, "Day10")

# Print the results for 2-day data
print(DAY10_results$summary)
print(DAY10_results$tukey)


## Fig S1E RPMI LP1 Knockdown lines
##------------------------------------------------------------------------------
# Define the function to perform ANOVA and Tukey HSD
perform_anova_tukey <- function(data, cell_type, treat_condition) {
  # Filter the data based on cell type and treatment condition
  filtered_data <- data %>%
    filter(Cell_type == cell_type & Treat == treat_condition)
  
  # Perform ANOVA
  res.aov <- aov(Cell_number ~ Genotype, data = filtered_data)
  aov_summary <- tidy(res.aov)
  
  # Perform Tukey HSD test
  tukey_results <- TukeyHSD(res.aov)
  
  list(summary = aov_summary, tukey = tukey_results)
}

# Load the data
RPMI_LP1_KD <- read_excel("Fig S1E RPMI LP1 KD.xlsx", sheet = "RPMI LP1 KD") %>%
  select(Cell_type, Genotype, Cell_number, Treat) %>%
  drop_na()

# Perform the analysis for RPMI_No_Dox treatment
RPMI_NoDox_results <- perform_anova_tukey(RPMI_LP1_KD, "RPMI8226", "No_Dox")
print(RPMI_NoDox_results$summary)
print(RPMI_NoDox_results$tukey)

# Perform the analysis for RPMI+Dox treatment
RPMI_Dox_results <- perform_anova_tukey(RPMI_LP1_KD, "RPMI8226", "Dox")
print(RPMI_Dox_results$summary)
print(RPMI_Dox_results$tukey)


# Perform the analysis for LP1_No_Dox treatment
LP1_NoDox_results <- perform_anova_tukey(RPMI_LP1_KD, "LP1", "No_Dox")
print(LP1_NoDox_results$summary)
print(LP1_NoDox_results$tukey)

# Perform the analysis for LP1+Dox treatment
LP1_Dox_results <- perform_anova_tukey(RPMI_LP1_KD, "LP1", "Dox")
print(LP1_Dox_results$summary)
print(LP1_Dox_results$tukey)


## Fig S8C FHA
#-------------------------------------------------------------------------------


FHA_lines <-  read_excel("Fig S8C FHA.xlsx", sheet = "020423") %>% 
  select(Genotype, Values, Day)

perform_anova_FHA <- function(data, day) {
  
  day_data <-  data %>% filter(Day == day)
  res.aov <- aov(Values ~ Genotype, data = day_data)
  
  res.aov_summary <- tidy(res.aov)
  tukey_results <- TukeyHSD(res.aov)
  
  list(summary = res.aov_summary, tukey = tukey_results)
}

#Perform anova and tukey on Day2

Day8 <- perform_anova_FHA(FHA_lines, "Day 8")

print(Day8$summary)
print(Day8$tukey)


###Fig S8 D Rescue_viability
#-------------------------------------------------------------------------------
viability <- read_excel("Fig S8 D Rescue_viability.xlsx")##Figure S8D
head(viability)

# Prepare data
df_viability <- viability %>%
  select(Genotypes, Count_number, `Total live cells`) %>%
  rename(Cell_count = `Total live cells`) %>%
  drop_na()

head(df_viability)

# Perform ANOVA
viability_aov <- aov(Cell_count ~ Genotypes, data = df_viability)

# Summary of the analysis
viability_aov_summary <- tidy(viability_aov)
viability_tukey <- TukeyHSD(viability_aov)

# Print ANOVA summary and Tukey HSD results
print(viability_aov_summary)
print(viability_tukey)

# Optional: Generate a report of the ANOVA (if the report package is installed)
report(viability_aov)


###RTPCR MAF targets

RT_PCR <- read_excel("..//..//RT_PCR/071924 maf targets raw results_v1.xls", sheet = "CCND2") 

names(RT_PCR)

RT_PCR <- RT_PCR %>% 
  select(Sample, Gene, Foldchange) %>% 
  drop_na()

head(RT_PCR)

# Perform one-way ANOVA
anova_result <- aov(Foldchange ~ Sample, data = RT_PCR)
summary(anova_result)

# Conduct Tukey's HSD test
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)

# Extract comparisons with `NON+DOX`
tukey_non_dox <- tukey_result$Sample[grep("NON\\+DOX", rownames(tukey_result$Sample)),]
tukey_non_dox

