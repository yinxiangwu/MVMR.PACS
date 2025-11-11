# Load necessary libraries
library(ggplot2)
library(reshape2)
library(readxl)
library(RColorBrewer)

# Read the correlation matrix from the Excel file
file_path <- "correlation_matrix.xlsx"  # Adjust if necessary
corr_matrix <- read_excel(file_path, col_names = TRUE)

# Convert to matrix format (remove the first column if it contains row names)
corr_matrix <- as.matrix(corr_matrix[,-1])  # Adjust if row names are in the first column
rownames(corr_matrix) <- colnames(corr_matrix)

# Convert the matrix into a long format for ggplot
corr_long <- melt(corr_matrix)

# Reverse the y-axis order to place Var1 in the upper left corner
corr_long$Var2 <- factor(corr_long$Var2, levels = rev(unique(corr_long$Var2)))

# Create the heatmap
jpeg('heat_map.jpeg',height = 4200, width = 4200, res = 600)
ggplot(data = corr_long, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white", size = 0.1) +  # Thin white gridlines
  scale_fill_gradientn(colors = brewer.pal(9, "RdBu"), limits = c(-1, 1), name = "Correlation") +  # Better color scale
  theme_classic() +
  labs(title = "Correlation Matrix for SNP-exposure associations", x = "", y = "") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.position = 'left',
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold")
  ) +
  coord_fixed()  # Ensure square aspect ratio
dev.off()


library(corrplot)
corrplot(corr_matrix, method = 'shade', title = 'Correlation Matrix for SNP-exposure associations', tl.pos = 'l')


# Load necessary libraries
library(ggplot2)
library(reshape2)
library(RColorBrewer)

# Create the data frame with the specified values
data_matrix <- data.frame(
  Var1 = "Column",  # Single column name
  Var2 = factor(1:10, levels = 10:1),  # Reverse order to place 1 at the top
  value = c(1, 1.2, 0.8, 0, 0, 0, 0, 0, 0.3, 0.6)
)

# Create the heatmap
jpeg('heat_map2.jpeg',height = 2800, width = 700, res = 600)
ggplot(data_matrix, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white", size = 0.1) +  # Thin white borders for clarity
  scale_fill_gradientn(
    colors = brewer.pal(9, "OrRd"),
    name = expression(bold(beta))  # Legend title as bold LaTeX \bm{\beta}
  ) +
  theme_void() +  # Removes all axes, backgrounds, and lines
  labs(title = "") +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    legend.position = "left"  # Keep the legend on the right for clarity
  ) +
  coord_fixed()  # Ensure square aspect ratio
dev.off()

# Create the data frame with the specified values
data_matrix <- data.frame(
  Var1 = "Column",  # Single column name
  Var2 = factor(1:10, levels = 10:1),  # Reverse order to place 1 at the top
  value = c(0.9, 0.9, 0.9, 0, 0, 0, 0, 0, 0.25, 0.5)
)

# Create the heatmap
jpeg('heat_map3.jpeg',height = 2800, width = 700, res = 600)
ggplot(data_matrix, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white", size = 0.1) +  # Thin white borders for clarity
  scale_fill_gradientn(
    colors = brewer.pal(9, "OrRd"),
    breaks = c(0,0.3,0.6,0.9)
  ) +
  theme_void() +  # Removes all axes, backgrounds, and lines
  labs(title = "") +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.position = "none"  # Keep the legend on the right for clarity
  ) +
  coord_fixed()  # Ensure square aspect ratio
dev.off()

