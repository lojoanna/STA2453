# 02_analysis_functions.R
# Custom functions for spatial transcriptomics analysis.
#
# This script defines custom functions used throughout the analysis,
# including coordinate rotation, tissue expression plotting, and GO enrichment.
#

#' Compute Pearson Correlations Between Matrix Columns
#'
#' If only `matrix1` is supplied, returns pairwise correlations among its columns.
#' If `matrix2` is supplied, returns correlations between columns of `matrix1` and `matrix2`.
#'
#' @param matrix1 A numeric matrix.
#' @param matrix2 An optional numeric matrix.
#'
#' @return A correlation matrix.
#' @export
correlation <- function(matrix1, matrix2 = NULL) {
  # Single-matrix correlations
  if (is.null(matrix2)) {
    output <- (t(scale(matrix1)) %*% scale(matrix1)) / (nrow(matrix1) - 1)
  } else {
    # Two-matrix correlations
    output <- (t(scale(matrix1)) %*% scale(matrix2)) / (nrow(matrix1) - 1)
  }
  
  # Return correlation matrix
  return(output)
}


#' Scale a Numeric Vector to [0, 1]
#'
#' Rescales `vector` so that its minimum maps to 0 and maximum maps to 1.
#'
#' @param vector A numeric vector.
#'
#' @return A numeric vector scaled between 0 and 1.
#' @export
scale_01 <- function(vector) {
  # Convert to plain numeric vector
  vec <- as.vector(vector)
  
  # Get minimum value
  min_value <- min(vec, na.rm = TRUE)
  
  # Scale to [0, 1]
  output <- (vec - min_value) / (max(vec, na.rm = TRUE) - min_value)
  
  # Return scaled vector
  return(output)
}


#' Rotate Coordinates for Visualization
#'
#' This function rotates x and y coordinates by a specified number of degrees,
#' with options to center, scale, and flip the coordinates. The function is used to restore
#' the original orientation of spatial transcriptomics data for proper visualization.
#'
#' @param x A numeric vector of x coordinates.
#' @param y A numeric vector of y coordinates.
#' @param n_degrees Number of degrees to rotate the coordinates.
#' @param center Logical; if TRUE, centers the coordinates.
#' @param scale Logical; if TRUE, scales the coordinates.
#' @param flip_x Logical; if TRUE, flips the x coordinates.
#' @param flip_y Logical; if TRUE, flips the y coordinates.
#'
#' @return A data frame with columns pos_x and pos_y representing the rotated coordinates.
#' @export
rotate_coordinates <- function(x, y, n_degrees = 0, center = FALSE, scale = FALSE, flip_x = FALSE, flip_y = FALSE) {
  # Convert degrees to radians
  theta <- n_degrees * (pi / 180)
  
  # Apply rotation matrix to x and y coordinates
  x_rotated <- (cos(theta) * x) - (sin(theta) * y)
  y_rotated <- (sin(theta) * x) + (cos(theta) * y)
  
  # Center and/or scale the coordinates
  x_rotated <- scale(x = x_rotated, center = center, scale = scale)
  y_rotated <- scale(x = y_rotated, center = center, scale = scale)
  
  # Flip x or y coordinates if required
  if (flip_x) {
    x_rotated <- -1 * x_rotated
  }
  if (flip_y) {
    y_rotated <- -1 * y_rotated
  }
  
  # Return output as a data frame
  output <- data.frame(pos_x = x_rotated, pos_y = y_rotated)
  return(output)
}


#' Cross-Expression Correlation Between Cell-Neighbor Pairs
#'
#' Identifies mutually exclusive expression pairs and computes Pearson
#' correlations between cells and their neighbors for each gene pair.
#'
#' @param data A cells-by-genes expression matrix.
#' @param locations A cells-by-coordinate(s) matrix or data frame.
#' @param neighbor Integer: the n-th nearest neighbor (1 = first).
#' @param output_matrix Logical: if TRUE, returns full correlation matrix.
#'
#' @return A data frame of gene pairs with correlation values, or (if
#'   `output_matrix = TRUE`) the full symmetric correlation matrix.
#' @export
cross_expression_correlation <- function(data, locations,
                                         neighbor = 1,
                                         output_matrix = FALSE) {
  # Compute neighbor indices
  neighbor_idx <- neighbor + 1
  nn <- RANN::nn2(locations, locations,
                  k = neighbor_idx,
                  searchtype = "standard",
                  treetype = "kd")
  distances <- nn$nn.idx[, neighbor_idx]
  data_temp <- data[distances, ]
  
  # Create mutually exclusive masks
  mask_data      <- (data > 0) * data
  mask_data_temp <- (data_temp > 0) * data_temp
  
  X <- mask_data * (1 - mask_data_temp)
  Y <- mask_data_temp * (1 - mask_data)
  
  # Keep original expression values where mask is TRUE
  X <- X * data
  Y <- Y * data_temp
  
  # Compute correlations and symmetrize
  corr_mat <- correlation(X, Y)
  corr_mat <- (corr_mat + t(corr_mat)) / 2
  
  if (output_matrix) {
    # Return full correlation matrix
    return(corr_mat)
  }
  
  # Return tidy edge list of upper triangle
  ids <- which(upper.tri(corr_mat), arr.ind = TRUE)
  result <- data.frame(
    gene1       = rownames(corr_mat)[ids[, 1]],
    gene2       = colnames(corr_mat)[ids[, 2]],
    correlation = upper_tri(corr_mat)
  )
  
  return(result)
}


#' Plot Tissue Expression Based on Two Genes
#'
#' This function creates a plot of cell locations colored based on the binary expression of two specified genes.
#' It can optionally highlight differences based on neighbor-specific expression, thereby revealing how the spatial context
#' of cells may influence gene expression. This is particularly useful for identifying patterns of coordinated expression
#' between adjacent cells.
#'
#' @param expr_data A matrix or data frame containing gene expression data. Columns should include the specified genes.
#' @param location_df A data frame with spatial coordinates; columns will be renamed to "x" and "y".
#' @param gene1 A character string specifying the first gene.
#' @param gene2 A character string specifying the second gene.
#' @param neighbor_expression Logical; if TRUE, highlights cells based on neighbor-specific expression.
#' @param neighbor Integer; the neighbor index to consider (0-based offset internally).
#' @param point_size Numeric; the size of the points in the plot.
#' @param scale_bar Numeric; length of the scale bar to annotate on the plot.
#'
#' @return A ggplot object representing the tissue expression plot.
#' @export
tissue_expression_plot <- function(expr_data, location_df, gene1, gene2, 
                                   neighbor_expression = TRUE, neighbor = 1, 
                                   point_size = 0, scale_bar = 0) {
  # Subset and binarize expression data for the two genes
  expr_data_sub <- expr_data[, c(gene1, gene2)]
  gene_names <- colnames(expr_data_sub)
  colnames(expr_data_sub) <- c("gene1", "gene2")
  expr_data_sub[expr_data_sub > 0] <- 1
  
  # Ensure location data has proper column names
  colnames(location_df) <- c("x", "y")
  
  # Compute scale bar coordinates based on quantiles
  x_start <- as.numeric(quantile(location_df$x, probs = 0.95))
  x_end   <- x_start + scale_bar
  y_start <- as.numeric(quantile(location_df$y, probs = 0))
  y_end   <- y_start
  
  # If not considering neighbor-specific expression, assign cell types based on expression
  if (!neighbor_expression) {
    cell_type <- vector(mode = "character", length = nrow(expr_data_sub))
    cell_type[expr_data_sub[, 1] == 1] <- "Gene1"
    cell_type[expr_data_sub[, 2] == 1] <- "Gene2"
    cell_type[(expr_data_sub[, 1] == 1) & (expr_data_sub[, 2] == 1)] <- "Both"
    cell_type[(expr_data_sub[, 1] == 0) & (expr_data_sub[, 2] == 0)] <- "Neither"
    
    location_df <- data.frame(location_df, type = cell_type)
    location_df$type <- factor(location_df$type, levels = c("Neither", "Both", "Gene2", "Gene1"))
    location_df <- location_df[order(location_df$type), ]
    
    plot_obj <- ggplot(location_df) +
      aes(x = x, y = y, color = type) +
      geom_point(size = point_size) +
      scale_color_manual(values = c("Neither" = "gray88", "Both" = "chartreuse3", 
                                    "Gene2" = "deepskyblue4", "Gene1" = "brown3"),
                         labels = c("Neither" = "Neither", "Both" = "Both", 
                                    "Gene2" = gene_names[2], "Gene1" = gene_names[1])) +
      labs(x = "X Coordinates", y = "Y Coordinates", color = "") +
      guides(color = guide_legend(override.aes = list(size = 2), reverse = TRUE)) +
      annotate("segment", x = x_start, xend = x_end, y = y_start, yend = y_end, 
               linewidth = 1, color = "black") +
      theme_classic() +
      theme(axis.text = element_text(size = 10),
            axis.title = element_text(size = 12),
            legend.text = element_text(size = 10))
    return(plot_obj)
  }
  
  # Else: consider neighbor-specific expression
  neighbor_index <- neighbor + 1
  distances <- RANN::nn2(location_df, location_df, k = neighbor_index,
                         searchtype = "standard", treetype = "kd")
  distances <- distances$nn.idx[, neighbor_index]
  neighbor_expr_df <- expr_data_sub[distances, ]
  
  cells <- rownames(expr_data_sub)
  neighbors <- rownames(expr_data_sub)[distances]
  
  # Determine which cell-neighbor pairs meet mutually exclusive expression criteria
  pair1 <- expr_data_sub[, 1] > 0 & expr_data_sub[, 2] == 0 &
    neighbor_expr_df[, 2] > 0 & neighbor_expr_df[, 1] == 0
  pair2 <- expr_data_sub[, 2] > 0 & expr_data_sub[, 1] == 0 &
    neighbor_expr_df[, 1] > 0 & neighbor_expr_df[, 2] == 0
  
  gene1_cells <- unique(c(cells[pair1], neighbors[pair2]))
  gene2_cells <- unique(c(cells[pair2], neighbors[pair1]))
  
  light <- vector(mode = "character", length = nrow(expr_data_sub))
  light[rownames(expr_data_sub) %in% gene1_cells] <- "Gene1"
  light[rownames(expr_data_sub) %in% gene2_cells] <- "Gene2"
  light[light == ""] <- "Neither"
  
  meta_df <- data.frame(location_df, light = light)
  meta_df$light <- factor(meta_df$light)
  
  meta_neither <- meta_df[meta_df$light == "Neither", ]
  meta_others  <- meta_df[meta_df$light != "Neither", ]
  meta_reordered <- rbind(meta_neither, meta_others)
  
  plot_obj <- ggplot(meta_reordered) +
    aes(x = x, y = y, color = light) +
    geom_point(size = point_size) +
    scale_color_manual(values = c("Neither" = "gray88", "Both" = "chartreuse3", 
                                  "Gene2" = "deepskyblue4", "Gene1" = "brown3"),
                       labels = c("Neither" = "Neither", "Both" = "Both", 
                                  "Gene2" = gene_names[2], "Gene1" = gene_names[1])) +
    labs(x = "X Coordinates", y = "Y Coordinates", color = "") +
    guides(color = guide_legend(override.aes = list(size = 2))) +
    annotate("segment", x = x_start, xend = x_end, y = y_start, yend = y_end, 
             linewidth = 1, color = "black") +
    theme_classic() +
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 10))
  return(plot_obj)
}


#' Compute GO Enrichment via Hypergeometric Test
#'
#' For a given set of test genes, performs GO term enrichment against a background
#' via hypergeometric tests, applies FDR correction, and returns significant terms.
#'
#' @param test_genes Character vector of gene symbols to test.
#' @param background_genes Character vector of background gene symbols; defaults to NULL (all genes).
#' @param species Character; "mouse" or "human". Default: "mouse".
#' @param GO_group_min_size Integer; minimum GO group size. Default: 20.
#' @param GO_group_max_size Integer; maximum GO group size. Default: 500.
#' @param alpha Numeric; FDR threshold for significance. Default: 0.05.
#'
#' @return A list with elements:
#'   - p_values: data frame of GO terms, p-values, and significance flags
#'   - GO_matrix: binary GO×gene matrix used for testing
#' @export
GO_enrichment <- function(test_genes, background_genes = NULL,
                          species = "mouse",
                          GO_group_min_size = 20,
                          GO_group_max_size = 500,
                          alpha = 0.05) {
  if (species == "mouse") {
    library(org.Mm.eg.db)
  }
  if (species == "human") {
    library(org.Hs.eg.db)
  }
  library(GO.db)
  conflicted::conflict_prefer("select", "dplyr")
  conflicted::conflict_prefer("filter", "dplyr")
  conflicted::conflict_prefer("rename", "dplyr")
  conflicted::conflict_prefer("as.matrix", "base")
  conflicted::conflict_prefer("union", "base")
  conflicted::conflict_prefer("intersect", "base")
  
  # Map GO IDs to gene IDs
  if (species == "mouse") {
    goAnnots <- as.list(org.Mm.egGO2ALLEGS)
  }
  if (species == "human") {
    goAnnots <- as.list(org.Hs.egGO2ALLEGS)
  }
  
  # Convert to long tibble of GO_id × NCBI_ID
  goAnnots <- tibble(GO_id = names(goAnnots), NCBI_ID = goAnnots) %>% unnest(NCBI_ID)
  
  # Map NCBI IDs to symbols
  if (species == "mouse") {
    goAnnots <- inner_join(goAnnots, org.Mm.egSYMBOL %>% as.data.frame() %>% rename(NCBI_ID = gene_id, gene_symbol = symbol))
  }
  if (species == "human") {
    goAnnots <- inner_join(goAnnots, org.Hs.egSYMBOL %>% as.data.frame() %>% rename(NCBI_ID = gene_id, gene_symbol = symbol))
  }
  
  # Remove duplicates
  goAnnots %<>% distinct()
  
  # Add in aspect
  goAnnots %<>% mutate(GO_aspect = Ontology(GO_id))
  
  # Add in GO group name
  goAnnots %<>% mutate(GO_name = Term(GO_id))
  
  # Add in GO group size info
  goAnnots %<>% inner_join(goAnnots %>% group_by(GO_id) %>% dplyr::count() %>% rename(GO_group_size = n))
  
  # Build GO × gene matrix
  go_wide <- goAnnots %>%
    select(GO_id, gene_symbol) %>%
    mutate(value = 1) %>%
    pivot_wider(
      names_from  = GO_id,
      id_cols     = gene_symbol,
      values_fill = 0
    )
  go_matrix <- as.matrix(go_wide %>% select(-gene_symbol))
  rownames(go_matrix) <- go_wide %>% pull(gene_symbol)
  
  # Convert annotations to data frame
  goAnnots <- as.data.frame(goAnnots)
  
  # Define background and test sets
  if (is.null(background_genes)) {
    background_genes <- rownames(go_matrix)
  }
  test_genes <- intersect(test_genes, rownames(go_matrix))
  background_genes <- intersect(background_genes, rownames(go_matrix))
  background_genes <- union(test_genes, background_genes)
  
  # Subset GO matrix to background genes and filter GO groups
  go_matrix <- go_matrix[rownames(go_matrix) %in% background_genes,]
  go_matrix <- go_matrix[, colSums(go_matrix) != 0]
  annot <- data.frame(GO_id = unique(goAnnots$GO_id), GO_name = unique(goAnnots$GO_name))
  annot <- annot[annot$GO_id %in% colnames(go_matrix),]
  colnames(go_matrix) <- annot$GO_name
  go_matrix <- t(go_matrix)
  
  # Curate GO groups
  GO_groups <- data.frame(group = rownames(go_matrix), size = as.numeric(rowSums(go_matrix)))
  
  if (is.null(GO_group_max_size)){GO_group_max_size = max(GO_groups$size)}
  
  go_matrix <- go_matrix[(GO_groups$size >= GO_group_min_size) & (GO_groups$size <= GO_group_max_size),]
  GO_groups <- GO_groups[(GO_groups$size >= GO_group_min_size) & (GO_groups$size <= GO_group_max_size),]
  
  # GO enrichment: hypergeometric test per GO group
  pvals <- vector(mode = "numeric", length = nrow(go_matrix))
  
  for (i in 1:nrow(go_matrix)){
    
    overlap  <- go_matrix[i,] == 1
    overlap  <- names(overlap[overlap])
    overlap  <- length(intersect(overlap, test_genes))
    
    k <- length(test_genes)    # Sample size: number of genes picked from GO group (balls picked)
    q <- overlap - 1           # Observed successes: number of genes from k actually observed in GO group (number of white balls seen)
    m <- sum(go_matrix[i,])    # Possible successes: number of all genes present in GO group (number of white balls in the population)
    n <- ncol(go_matrix) - m   # Possible failures:  number of all genes absent in GO group (number of black balls in the population)
    
    pvals[i] <- phyper(q = q, m = m, n = n, k = k, lower.tail = FALSE)
  }
  
  # Process p-values
  pvals_orig <- pvals
  pvals_FDR <- p.adjust(pvals, method = "BH")
  pvals <- data.frame(GO_groups, pvals_orig, pvals_FDR)
  
  pvals_orig[pvals_orig == 0] <- min(pvals_orig[pvals_orig > 0])
  pvals_FDR[pvals_FDR == 0]   <- min(pvals_FDR[pvals_FDR > 0])
  
  pvals <- data.frame(pvals, sig_GO <- as.numeric(pvals$pvals_FDR <= alpha))
  
  # Append full GO annotations
  pvals <- data.frame(GO_id = goAnnots[match(pvals$group, goAnnots$GO_name), 1], 
                      group = pvals$group, 
                      GO_aspect = goAnnots[match(pvals$group, goAnnots$GO_name), 4], 
                      pvals[,2:ncol(pvals)])
  rownames(pvals) <- 1:nrow(pvals)
  rownames(go_matrix) <- pvals$GO_id
  
  # Combine output
  output <- list(pvals, go_matrix)
  names(output) <- c("p_values","GO_matrix")
  
  return(output)
}