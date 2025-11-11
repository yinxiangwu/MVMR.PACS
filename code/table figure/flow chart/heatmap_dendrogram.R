set.seed(123)

num_runs <- 100
num_vars <- 32
variables <- traits

enough_iv_strength <- 
  
  # Simulate clustering results
  clustering_results <- t(sapply(res_repeat_v2$grouping_hist,function(x) {as.numeric(strsplit(x, "-")[[1]])}))
rownames(clustering_results) <- NULL

n_cls <- apply(clustering_results,1,function(x) length(unique(x)))

# Convert to data frame for convenience
clustering_results_df <- as.data.frame(clustering_results)
colnames(clustering_results_df) <- paste0("Var", 1:num_vars)

# Simulate binary significance results (1 significant, 0 not significant)
significance_df <- ifelse(res_repeat_v2$beta.hist == 0, 0, (res_repeat_v2$beta.hist)*(2*pnorm(abs(res_repeat_v2$beta.hist/res_repeat_v2$se.hist),lower.tail = FALSE) < 0.05))
colnames(significance_df) <- variables
rownames(significance_df) <- paste0("Run", 1:num_runs)

library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

# Co-occurrence matrix
co_occurrence <- matrix(0, num_vars, num_vars)
colnames(co_occurrence) <- rownames(co_occurrence) <- variables

for(i in 1:num_vars) {
  for(j in 1:num_vars) {
    co_occurrence[i, j] <- mean(clustering_results[, i] == clustering_results[, j])
  }
}

# Hierarchical clustering from co-occurrence
dist_matrix <- as.dist(1 - co_occurrence)
hc <- hclust(dist_matrix, method="average")

# (3) Plot heatmap preserving original order with dendrogram
sig_matrix <- as.matrix(significance_df)

# sig_matrix <- sig_matrix[res_repeat_v2$iv.hist.dt1 > 7, ]

# Important: Define dendrogram separately without affecting column order
column_dend <- as.dendrogram(hc)

# Define color scale: blue = negative, white = 0, red = positive
col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

# Plot the heatmap
jpeg('~/OneDrive - UW/UW biostats/2023-24 RA/MR with highly correlated exposure/manuscript/table figure/heatmap_5e08_0830.jpeg',width = 3200,height = 2400,res = 600)
Heatmap(sig_matrix,
        name = "Beta",
        col = col_fun,
        show_row_names = FALSE,
        cluster_rows = TRUE,
        show_row_dend = FALSE,
        cluster_columns = column_dend,
        column_order = variables,
        column_dend_height = unit(3, "cm"),
        column_names_side = "top",
        heatmap_legend_param = list(title = "beta"))
dev.off()