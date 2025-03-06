# 02_analysis_functions.R
# Custom functions for spatial transcriptomics analysis.
#
# This script defines custom functions used throughout the analysis,
# including coordinate rotation, tissue expression plotting, and co-occurrence statistics.
#

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

#' Compute Co-occurrence Statistics Between Two Matrices
#'
#' This function computes the number of shared elements between the columns of two matrices.
#' It optionally converts the matrices to a sparse format and binarizes them (all positive values are set to 1).
#'
#' @param matrix1 A numeric matrix.
#' @param matrix2 A numeric matrix.
#' @param sparse Logical; if TRUE, converts matrices to sparse format.
#' @param binarize Logical; if TRUE, binarizes the matrices (all positive values are set to 1).
#'
#' @return A list with two elements: 'cooccur' (the co-occurrence counts) and 'cells'
#'   (the dot product of matrix1 with itself).
#' @export
get_cooccurrence_stats <- function(matrix1, matrix2, sparse = TRUE, binarize = TRUE) {
  # Convert to sparse matrices if specified
  if (sparse) {
    matrix1 <- Matrix(data = Matrix::as.matrix(matrix1), sparse = TRUE)
    matrix2 <- Matrix(data = Matrix::as.matrix(matrix2), sparse = TRUE)
  }
  
  # Binarize matrices if requested (all positive values become 1)
  if (binarize) {
    matrix1[matrix1 > 0] <- 1
    matrix2[matrix2 > 0] <- 1
  }
  
  # Compute the number of shared elements between columns
  cooccur <- t(matrix1) %*% matrix2
  cells   <- t(matrix1) %*% matrix1
  
  output <- list(cooccur = cooccur, cells = cells)
  return(output)
}