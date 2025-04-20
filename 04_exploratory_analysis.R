# 04_exploratory_analysis.R
# Exploratory data analysis for spatial transcriptomics data.
#
# This script performs exploratory analyses, including:
# - Visualization of mRNA count distributions and cell volume distributions.
# - Analysis of the relationship between total mRNA counts and cell volume.
# - Gene-gene co-expression analysis across slices.
# - Field-of-view (FOV) analysis and nearest-neighbor relationships.
# - Visualization of cell-neighbor gene expression similarity.
#
# It assumes that the data has been preprocessed in 03_preprocessing.R and that custom functions are available.
#

# Plot distribution of total mRNA counts per cell
counts_df <- data.frame()
for (i in seq_along(data_list)) {
  total_counts_df <- data.frame(counts = rowSums(data_list[[i]]),
                                sample = names(data_list)[i])
  counts_df <- rbind(counts_df, total_counts_df)
}

ggplot(counts_df) +
  aes(x = counts, y = sample, fill = sample) +
  geom_density_ridges(scale = 1, alpha = 0.8) +
  scale_x_log10() +
  geom_vline(xintercept = 55, linetype = "dotted", color = "brown2") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Total Number of mRNA Counts per Cell", y = "",
       subtitle = "Dotted Line = Quality Control Threshold (55)") +
  scale_fill_viridis_d()

# Plot distribution of cell volume per slice
cell_volume_df <- data.frame()
for (i in seq_along(metadata_list)) {
  temp_volume_df <- data.frame(volume = metadata_list[[i]]$volume,
                               sample = names(data_list)[i])
  cell_volume_df <- rbind(cell_volume_df, temp_volume_df)
}

ggplot(cell_volume_df) +
  aes(x = volume, y = sample, fill = sample) +
  geom_density_ridges(scale = 1, alpha = 0.8) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Cell Volume", y = "") +
  scale_fill_viridis_d()

# Analyze total mRNA counts versus cell volume with heteroscedasticity test
counts_vs_volume_df <- data.frame()
for (i in seq_along(metadata_list)) {
  temp_df <- data.frame(counts = rowSums(data_list[[i]]),
                        volume = metadata_list[[i]]$volume)
  counts_vs_volume_df <- rbind(counts_vs_volume_df, temp_df)
}

# Perform Breusch-Pagan test for heteroscedasticity
test <- bptest(lm(counts_vs_volume_df$counts ~ counts_vs_volume_df$volume))

ggplot(counts_vs_volume_df) +
  aes(x = volume, y = counts) +
  geom_point(size = 0, alpha = 0.05) +
  labs(x = "Cell Volume",
       y = "Total Number of mRNA Counts per Cell",
       subtitle = str_c("Spearman's Correlation = ", 
                        signif(cor(counts_vs_volume_df, method = "spearman")[1, 2], digits = 2),
                        "\nBreusch-Pagan Test for Homoscedasticity, p-value ≤ 2.2e-16")) +
  theme_classic()

# Gene-gene co-expression analysis between slices
co_expression <- matrix(data = 0,
                        nrow = choose(ncol(data_list[[1]]), 2),
                        ncol = length(data_list))
colnames(co_expression) <- names(data_list)

for (i in seq_along(data_list)) {
  co_expression[, i] <- upper_tri(cor(as.matrix(data_list[[i]])))
  print(i / length(data_list))
}

# Plot gene-gene co-expression example: one slice versus another
temp_co_expr <- co_expression[, c(5, 8)]
colnames(temp_co_expr) <- c("x", "y")

ggplot(temp_co_expr) +
  aes(x = x, y = y) +
  geom_point(size = 0, alpha = 0.8) +
  theme_classic() +
  geom_abline(intercept = 0, slope = 1, color = "brown3", linewidth = 0.5) +
  labs(x = "Gene-Gene Co-expression\nBrain 2, Middle Slice",
       y = "Gene-Gene Co-expression\nBrain 3, Middle Slice",
       subtitle = str_c("Spearman's Correlation = ", 
                        signif(cor(temp_co_expr, method = "spearman")[1, 2], digits = 2)))

# Generate heatmap of gene-gene co-expression correlations
co_expr_corr <- cor(co_expression, method = "spearman")
pheatmap(co_expr_corr)

# Generate gene-gene co-expression heatmap with hierarchical clustering to reveal modules
gene_corr_matrix <- cor(as.matrix(data_list[[5]]))
Heatmap(gene_corr_matrix, 
        row_names_gp = gpar(fontsize = 1),
        column_names_gp = gpar(fontsize = 1),
        column_title = "Genes", column_title_side = "bottom",
        row_title = "Genes", row_title_side = "right")

# Field-of-view (FOV) analysis and nearest-neighbor relationships
temp_metadata <- metadata_list[[2]]
fov_sample <- sample(x = unique(temp_metadata$fov), 
                     size = round(length(unique(temp_metadata$fov)) * 0.05))
temp_metadata$show_fov <- as.numeric(temp_metadata$fov %in% fov_sample)

ggplot(temp_metadata) +
  aes(x = x, y = y, color = show_fov) +
  geom_point(size = 0) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "", y = "", subtitle = "5% of FOVs Shown")

ggplot(temp_metadata) +
  aes(x = x, y = y, color = fov) +
  geom_point(size = 0) +
  theme_classic() +
  labs(x = "", y = "", color = "FOV ID") +
  theme(legend.position = "none")

# Count number of cells per FOV
fov_counts <- c()
for (i in seq_along(metadata_list)) {
  fov_counts <- c(fov_counts, as.numeric(table(metadata_list[[i]]$fov)))
}
fov_counts <- data.frame(fov = fov_counts)

ggplot(fov_counts) +
  aes(x = fov, y = after_stat(scaled)) +
  geom_density(fill = "lightblue", color = "darkblue") +
  labs(x = "Number of Cells within a Field-of-View (FOV)", y = "Density") +
  theme_classic()

# Nearest-neighbor FOV relationship analysis
neighbor_prop <- c()
for (i in seq_along(metadata_list)) {
  temp_metadata <- metadata_list[[i]]
  neighbor_metadata <- temp_metadata[
    RANN::nn2(data = temp_metadata[, c("x", "y")],
              query = temp_metadata[, c("x", "y")],
              k = 2,
              treetype = "kd",
              searchtype = "standard")$nn.idx[, 2],
  ]
  neighbor_prop[i] <- sum(temp_metadata$fov == neighbor_metadata$fov) / nrow(temp_metadata)
}
print(str_c("Proportion of Cells Whose Neighbors Are in Different FOV = ", 
            signif(1 - mean(neighbor_prop), digits = 2)))

# Analyze cell-neighbor gene expression similarity
cell_neighbor_corr <- c()
for (i in seq_along(metadata_list)) {
  temp_metadata <- metadata_list[[i]]
  neighbor_metadata <- temp_metadata[
    RANN::nn2(data = temp_metadata[, c("x", "y")],
              query = temp_metadata[, c("x", "y")],
              k = 2,
              treetype = "kd",
              searchtype = "standard")$nn.idx[, 2],
  ]
  
  colnames(neighbor_metadata)[colnames(neighbor_metadata) == "V1"] <- "cell_id"
  
  temp_expr <- data_list[[i]]
  neighbor_expr <- data_list[[i]][match(neighbor_metadata$cell_id, rownames(data_list[[i]])), , drop = FALSE]
  
  # Z-score normalize the expression data for correlation calculation
  temp_expr <- scale(t(temp_expr))
  neighbor_expr <- scale(t(neighbor_expr))
  
  output_corr <- as.numeric(colSums(temp_expr * neighbor_expr) / (nrow(temp_expr) - 1))
  cell_neighbor_corr <- c(cell_neighbor_corr, output_corr)
  
  print(i / length(metadata_list))
}
cell_neighbor_corr_df <- data.frame(neighbor_corr = cell_neighbor_corr)
ggplot(cell_neighbor_corr_df) +
  aes(x = neighbor_corr, y = after_stat(scaled)) +
  geom_density(fill = "lightblue", color = "darkblue") +
  labs(x = "Correlation Between Cells' and Their Neighbors' Gene Expression", y = "Density") +
  theme_classic()

# Example usage of tissue expression plotting function
# (Assumes that the custom function is defined in 02_analysis_functions.R)
tissue_expression_plot(expr_data = data_list[[2]],
                       location_df = metadata_list[[2]][, c("x", "y")],
                       gene1 = "Lgr6", gene2 = "Adra2b",
                       neighbor_expression = FALSE)
tissue_expression_plot(expr_data = data_list[[2]],
                       location_df = metadata_list[[2]][, c("x", "y")],
                       gene1 = "Lgr6", gene2 = "Adra2b",
                       neighbor_expression = TRUE)