# 03_preprocessing.R
# Preprocess spatial transcriptomics data.
#
# This script performs data cleaning and transformation steps, including:
# - Rotating spatial coordinates to restore original slice orientation.
# - Filtering low-quality cells based on mRNA counts.
# - Computing average gene expression across cells.
#
# It assumes that 'data_list' and 'metadata_list' have been loaded from 01_data_loading.R.
# It also assumes that custom functions are available from 03_analysis_functions.R.
#
# Source custom analysis functions
source("~/Documents/***COURSEWORK***/STA2453H/Research/(+) GitHub/scripts/02_analysis_functions.R")

# Define rotation degrees for counter-clockwise rotation
rotation_degrees <- c(175, 165, 234, 238, 41, 165, 224, 353, 65)
#rotation_degrees <- 360 - rotation_degrees

# Rotate each slice and display before and after plots
for (i in seq_along(metadata_list)) {
  temp_metadata <- metadata_list[[i]]
  
  # Plot original coordinates
  plot_obj <- ggplot(temp_metadata) +
    aes(x = x, y = y) +
    geom_point(size = 0) +
    labs(x = "", y = "") +
    theme_classic()
  print(plot_obj)
  
  # Apply rotation to coordinates using custom function
  rotated_coords <- rotate_coordinates(x = temp_metadata$x, y = temp_metadata$y,
                                       n_degrees = rotation_degrees[i])
  temp_metadata$x <- rotated_coords$pos_x
  temp_metadata$y <- rotated_coords$pos_y
  
  # Plot rotated coordinates
  plot_obj <- ggplot(temp_metadata) +
    aes(x = x, y = y) +
    geom_point(size = 0) +
    labs(x = "", y = "") +
    theme_classic()
  print(plot_obj)
  
  # Update metadata_list with rotated coordinates
  metadata_list[[i]] <- temp_metadata
}

# Filter cells and compute average gene expression per sample
avg_expression <- matrix(data = 0,
                         nrow = ncol(data_list[[1]]),
                         ncol = length(data_list),
                         dimnames = list(colnames(data_list[[1]]), names(data_list)))

for (i in seq_along(data_list)) {
  temp_expr <- data_list[[i]]
  temp_metadata <- metadata_list[[i]]
  
  # Remove cells with fewer than 50 counts (quality control)
  hits <- rowSums(temp_expr) >= 50
  temp_expr <- temp_expr[hits, ]
  temp_metadata <- temp_metadata[hits, ]
  
  # Compute average counts per gene
  avg_expression[, i] <- as.numeric(colMeans(temp_expr))
  
  # Update stored data
  data_list[[i]] <- temp_expr
  metadata_list[[i]] <- temp_metadata
}