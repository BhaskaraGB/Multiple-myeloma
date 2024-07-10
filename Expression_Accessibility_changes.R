# Load necessary packages
library(ggplot2)

# Prepare the data for all four sets (3 DAY)
data <- data.frame(
  Category = c("Decreased Accessibility Only", "Overlap", "Upregulated Genes Only", 
               "Decreased Accessibility Only", "Overlap", "Downregulated Genes Only", 
               "Increased Accessibility Only", "Overlap", "Downregulated Genes Only",
               "Increased Accessibility Only", "Overlap", "Upregulated Genes Only"),
  Count = c(7734, 210, 1021, 7460, 484, 816, 4914, 228, 1072, 4809, 333, 898),
  Set = c("Upregulated Genes", "Upregulated Genes", "Upregulated Genes", 
          "Decreased Accessibility with Downregulated Genes", "Decreased Accessibility with Downregulated Genes", "Decreased Accessibility with Downregulated Genes", 
          "Increased Accessibility with Downregulated Genes", "Increased Accessibility with Downregulated Genes", "Increased Accessibility with Downregulated Genes",
          "Increased Accessibility with Upregulated Genes", "Increased Accessibility with Upregulated Genes", "Increased Accessibility with Upregulated Genes")
)

# Convert 'Set' to factor to maintain order
data$Set <- factor(data$Set, levels = c("Upregulated Genes", "Decreased Accessibility with Downregulated Genes", "Increased Accessibility with Downregulated Genes", "Increased Accessibility with Upregulated Genes"))

# Ensure the levels of Category for proper stacking
data$Category <- factor(data$Category, levels = c("Overlap", "Upregulated Genes Only",
                                                  "Downregulated Genes Only","Decreased Accessibility Only","Increased Accessibility Only"))

# Load necessary packages
library(ggplot2)

cols <- c("Overlap" = "purple", "Upregulated Genes Only" = "red",
          "Downregulated Genes Only" = "blue","Decreased Accessibility Only" = "lightskyblue","Increased Accessibility Only" = "pink")


# Create the stacked bar chart with reduced gap
ATAC_expression_3d <- ggplot(data, aes(x = Set, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 1, color = "black") +
  geom_text(aes(label = Count), position = position_stack(vjust = 0.5), size = 6, fontface = "bold", color = "white") +
  scale_y_continuous(expand = expansion(0), limits = c(0, 11000)) +
  scale_fill_manual(values = cols) +
  theme_minimal() +
  labs(title = "Differential Gene Expression and Accessibility Genes (3d)", x = NULL, y = "Gene count") +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    legend.key.width = unit(1.5, "line"),
    legend.key.height = unit(1.5, "line"),
    legend.key.spacing.y = unit(1, "line"),
    legend.text = element_text(size = 16))

# Display the plot
print(ATAC_expression_3d)

# Save the plot
directory <- "/Users/bhaskara/Documents/Candy/"
filename <- paste0(directory, "/ATAC_expression_3d.tiff")

ggsave(filename, plot = ATAC_expression_3d, width = 10, height = 8, dpi = 300)


#___________________________________________________________________________________________________________________________________________________________

# Prepare the data for all four sets (5 DAY)
data <- data.frame(
  Category = c("Decreased Accessibility Only", "Overlap", "Upregulated Genes Only", 
               "Decreased Accessibility Only", "Overlap", "Downregulated Genes Only", 
               "Increased Accessibility Only", "Overlap", "Downregulated Genes Only",
               "Increased Accessibility Only", "Overlap", "Upregulated Genes Only"),
  Count = c(8568, 391, 1178, 8294, 665, 1023, 6574, 320, 1368, 6353, 541, 1028),
  Set = c("Upregulated Genes", "Upregulated Genes", "Upregulated Genes", 
          "Decreased Accessibility with Downregulated Genes", "Decreased Accessibility with Downregulated Genes", "Decreased Accessibility with Downregulated Genes", 
          "Increased Accessibility with Downregulated Genes", "Increased Accessibility with Downregulated Genes", "Increased Accessibility with Downregulated Genes",
          "Increased Accessibility with Upregulated Genes", "Increased Accessibility with Upregulated Genes", "Increased Accessibility with Upregulated Genes")
)

# Convert 'Set' to factor to maintain order
data$Set <- factor(data$Set, levels = c("Upregulated Genes", "Decreased Accessibility with Downregulated Genes", "Increased Accessibility with Downregulated Genes", "Increased Accessibility with Upregulated Genes"))

# Ensure the levels of Category for proper stacking
data$Category <- factor(data$Category, levels = c("Overlap", "Upregulated Genes Only",
                                                  "Downregulated Genes Only","Decreased Accessibility Only","Increased Accessibility Only"))

# Define colors
cols <- c("Overlap" = "purple", "Upregulated Genes Only" = "red",
          "Downregulated Genes Only" = "blue","Decreased Accessibility Only" = "lightskyblue","Increased Accessibility Only" = "pink")

# Create the stacked bar chart with reduced gap
ATAC_expression_5d <- ggplot(data, aes(x = Set, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 1, color = "black") +
  geom_text(aes(label = Count), position = position_stack(vjust = 0.5), size = 6, fontface = "bold", color = "white") +
  scale_y_continuous(expand = expansion(0), limits = c(0, 11000)) +
  scale_fill_manual(values = cols) +
  theme_minimal() +
  labs(title = "Differential Gene Expression and Accessibility Genes (5d)", x = NULL, y = "Gene count") +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    legend.key.width = unit(1.5, "line"),
    legend.key.height = unit(1.5, "line"),
    legend.key.spacing.y = unit(1, "line"),
    legend.text = element_text(size = 16)
  )

# Display the plot
print(ATAC_expression_5d)

# Save the plot
directory <- "/Users/bhaskara/Documents/Candy/"
filename <- paste0(directory, "/ATAC_expression_5d.tiff")

ggsave(filename, plot = ATAC_expression_5d, width = 10, height = 8, dpi = 300)


