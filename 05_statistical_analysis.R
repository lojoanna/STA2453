# 05_statistical_analysis.R
# Perform statistical and network analyses on spatial transcriptomics data.
#
# This script computes:
#  - Average gene expression correlations across slices
#  - Gene–gene co‐expression across slices and heatmap visualization
#  - Cell–neighbor correlation networks and heatmaps
#  - Meta‐analytic network consensus clustering and GO enrichment
#  - ROC‐based evaluation of correlation vs. p‐value predictions
#
# Assumes `data_list` and `metadata_list` are created by 01–03, and that
# custom functions are available from 02_analysis_functions.R.
#

# Source custom analysis functions
source("~/Documents/***COURSEWORK***/STA2453H/Research/(+) GitHub/scripts/02_analysis_functions.R")


# 1. Average gene expression comparisons ---------------------------------------

# Compute per‐slice average expression (cells with ≥55 counts)
avg_exp <- matrix(
  0,
  nrow = ncol(data_list[[1]]),
  ncol = length(data_list),
  dimnames = list(colnames(data_list[[1]]), names(data_list))
)

for (i in seq_along(data_list)) {
  temp <- data_list[[i]]
  keep <- rowSums(temp) >= 50
  temp <- temp[keep, ]
  avg_exp[, i] <- colMeans(temp)
}

# Example scatter: slice 1 vs slice 4
example_avg <- avg_exp[, c(1, 4)]
colnames(example_avg) <- c("x", "y")

ggplot(as.data.frame(example_avg)) +
  aes(x = x, y = y) +
  geom_point(size = 1, alpha = 0.8) +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(slope = 1, intercept = 0, color = "brown3", linewidth = 0.5) +
  theme_classic() +
  labs(
    x = "Average gene counts\nBrain 1, posterior slice",
    y = "Average gene counts\nBrain 2, posterior slice",
    subtitle = str_c(
      "Spearman's correlation = ",
      signif(cor(example_avg, method = "spearman")[1, 2], 2)
    )
  )

# Heatmap of average‐expression correlations
avg_corr <- cor(avg_exp, method = "spearman")
Heatmap(
  avg_corr,
  heatmap_legend_param = list(title = "Spearman ρ"),
  column_title = "Average Expression Correlations"
)


# 2. Gene–gene co‐expression across slices -------------------------------------

plan(multisession, workers = parallel::detectCores() - 1)

co_expr <- do.call(
  cbind,
  future_lapply(
    data_list,
    function(x) upper_tri(cor(as.matrix(x))),
    future.seed = TRUE
  )
)
colnames(co_expr) <- names(data_list)

# Example scatter: co‐expression slice 5 vs slice 8
example_co <- co_expr[, c(5, 8)]
colnames(example_co) <- c("x", "y")

ggplot(as.data.frame(example_co)) +
  aes(x = x, y = y) +
  geom_point(size = 0, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, color = "brown3", linewidth = 0.5) +
  theme_classic() +
  labs(
    x = "Co‐expression\nBrain 2, middle slice",
    y = "Co‐expression\nBrain 3, middle slice",
    subtitle = str_c(
      "Spearman's correlation = ",
      signif(cor(example_co, method = "spearman")[1, 2], 2)
    )
  )

# Heatmap of co‐expression correlations
co_corr <- cor(co_expr, method = "spearman")
Heatmap(
  co_corr,
  heatmap_legend_param = list(title = "Spearman ρ"),
  column_title = "Co‐expression Correlations"
)


# 3. Cell–neighbor correlation networks ---------------------------------------

plan(multisession, workers = parallel::detectCores() - 1)

options(future.globals.maxSize = 10 * 1024^3)

net_list <- future_lapply(
  seq_along(data_list),
  function(i) {
    cross_expression_correlation(
      data        = data_list[[i]],
      locations   = metadata_list[[i]][, c("x", "y")],
      output_matrix = TRUE
    )
  },
  future.seed = TRUE
)
network <- sapply(net_list, upper_tri)
colnames(network) <- names(data_list)


# Example scatter: neighbor correlation slice 1 vs slice 4
example_net <- network[, c(1, 4)]
colnames(example_net) <- c("x", "y")

ggplot(as.data.frame(example_net)) +
  aes(x = x, y = y) +
  geom_point(size = 1, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, color = "brown3", linewidth = 0.5) +
  theme_classic() +
  labs(
    x = "Neighbor correlation\nBrain 1, posterior slice",
    y = "Neighbor correlation\nBrain 2, posterior slice",
    subtitle = str_c(
      "Spearman's correlation = ",
      signif(cor(example_net, method = "spearman")[1, 2], 2)
    )
  )

# Heatmap of neighbor‐correlation correlations
net_corr <- cor(network, method = "spearman")
Heatmap(
  net_corr,
  heatmap_legend_param = list(title = "Spearman ρ"),
  column_title = "Cell–Neighbor Correlation Correlations"
)


# 4. Meta‐analytic network & community detection ------------------------------

# Build consensus matrix from upper‐triangular of `network`
med_vals <- matrixStats::rowMedians(network)

n_genes <- ncol(data_list[[1]])
meta_mat <- matrix(
  0,
  nrow = n_genes, ncol = n_genes,
  dimnames = list(colnames(data_list[[1]]), colnames(data_list[[1]]))
)
meta_mat[upper.tri(meta_mat)] <- med_vals
meta_mat <- (meta_mat + t(meta_mat)) / 2
diag(meta_mat) <- 0

# Convert to igraph
edges <- which(upper.tri(meta_mat), arr.ind = TRUE)
edge_df <- data.frame(
  gene1 = rownames(meta_mat)[edges[, 1]],
  gene2 = colnames(meta_mat)[edges[, 2]],
  weight = meta_mat[upper.tri(meta_mat)]
)
g <- graph_from_data_frame(edge_df, directed = FALSE)

# Leiden clustering, consensus over `iter` runs
iter <- 100
partitions <- replicate(
  iter,
  membership(cluster_leiden(
    g,
    weights = abs(E(g)$weight),
    objective_function = "modularity",
    n_iterations = 100
  ))
)

# Compute ARI matrix to pick best consensus partition
ARI <- matrix(0, nrow = iter, ncol = iter)
for (i in seq_len(iter)) {
  for (j in seq_len(iter)) {
    ARI[i, j] <- adjustedRandIndex(partitions[, i], partitions[, j])
  }
}

best_run <- which.max(matrixStats::colMedians(ARI))
clusters <- partitions[, best_run]

# Assign GO term labels to each cluster based on the most significantly enriched GO group
unique_clusts <- sort(unique(clusters))
go_labels <- data.frame()
for (k in 1:max(clusters)) {
  res <- GO_enrichment(
    test_genes       = names(clusters)[clusters == k],
    background_genes = names(clusters)
  )
  res <- res$p_values
  top_term <- res[as.logical(res$sig_GO),]
  if (nrow(top_term) == 0){print(k/max(clusters)); next}
  go_labels <- rbind(go_labels, data.frame(GO_group = top_term$group, cluster = k))
  print(k/max(clusters))
}

# Merge clusters and provide labels
label_clusters <- data.frame(genes = names(clusters), comm = as.numeric(clusters))
label_clusters[label_clusters$comm == 1 | label_clusters$comm == 2, ]$comm <- "Cell-cell signaling by brain cells"
label_clusters[label_clusters$comm == 4 | label_clusters$comm == 3, ]$comm <- "Cell identity/type markers"
label_clusters[label_clusters$comm == 5, ]$comm <- "Cell-cell signaling by immune cells"

# Annotate heatmap by community
group_vec    <- setNames(label_clusters$comm, label_clusters$genes)
group_levels <- sort(unique(group_vec))
group_colors <- setNames(c("chartreuse3", "deepskyblue4", "brown3"), group_levels)

row_anno <- HeatmapAnnotation(
  Community = as.factor(group_vec),
  col = list(Community = group_colors),
  which = "row",
  show_legend = TRUE,
  show_annotation_name = FALSE,
  annotation_legend_param = list(title = "Biological Community"))

col_anno <- HeatmapAnnotation(
  Community = as.factor(group_vec),
  col = list(Community = group_colors),
  which = "column",
  show_legend = FALSE,
  show_annotation_name = FALSE)

Heatmap(meta_mat,
        row_names_gp = gpar(fontsize = 1),
        column_names_gp = gpar(fontsize = 1),
        heatmap_legend_param = list(title = "Correlation"),
        top_annotation = col_anno,
        left_annotation = row_anno,
        column_title = "Genes",
        column_title_side = "bottom",
        row_title = "Genes",
        row_title_side = "right")


# 5. ROC evaluation of correlation → p‐value -------------------------------

# Pre‐compute observed correlations
obs_corr <- network

# Compute p‐values via permutation (e.g. 1000 shuffles)
iter_perm <- 1000
pvals <- matrix(
  0,
  nrow = nrow(obs_corr),
  ncol = ncol(obs_corr),
  dimnames = dimnames(obs_corr)
  )

# Original test statistics
test_stat <- obs_corr

# Permute cell labels and recompute neighbor correlations
for (i in seq_len(length(data_list))) {
  for (k in seq_len(iter_perm)) {
    perm_data <- data_list[[i]][sample(nrow(data_list[[i]])), ]
    tmp_corr <- cross_expression_correlation(
      perm_data,
      metadata_list[[i]][, c("x", "y")],
      output_matrix = FALSE
    )$correlation
    pvals[, i] <- pvals[, i] + (test_stat[, i] >= tmp_corr)
  }
  pvals[, i] <- 1 - (pvals[, i] / iter_perm)
  pvals[, i] <- p.adjust(pvals[, i], method = "BH")
}

pvals <- (pvals <= 0.05) * 1

# Compute TPR and FPR using correlations and p-values
outcome <- data.frame()

for (i in 1:ncol(test_stat)) {
  for (j in 1:ncol(pvals)) {
    x <- pROC::roc(response = pvals[, i], predictor = test_stat[, j])
    x <- data.frame(TPR = x$sensitivities, FPR = 1 - x$specificities,
                    data1 = colnames(pvals)[i], data2 = colnames(test_stat)[j])
    outcome <- rbind(outcome, x)
    print(str_c("i = ", i, "; j = ", j))
  }
}

# Assign labels (self, homologous, non-homologous) based on slice relationships
# and sample subset (100,000 points) of `outcome` so data can be plotted
result <- outcome
result <- data.frame(result, type = "X")
result <- result[sample(1:nrow(result), size = 100000, replace = FALSE), ]
result$type <- ifelse(
  result$data1 == result$data2, "Self", 
  ifelse((grepl("anterior",  result$data1) & grepl("anterior",  result$data2)) | 
         (grepl("middle",    result$data1) & grepl("middle",    result$data2)) | 
         (grepl("posterior", result$data1) & grepl("posterior", result$data2)),
         "Homologous",
         "Non-homologous"
  )
)
result$type <- factor(result$type, levels = rev(c("Self", "Homologous", "Non-homologous")))

# Plot AUROC
ggplot(result) + 
  aes(x = FPR, y = TPR, color = type) + 
  geom_point(size = 0) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") + 
  theme_minimal() +
  labs(
    x = "False positive rate (FPR)",
    y = "True positive rate (TPR)",
    color = "Slice relationships",
    subtitle = "ROC using correlations to predict p-values") +
  guides(color = guide_legend(reverse = TRUE, override.aes = list(size = 3))) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 10),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12)
  )