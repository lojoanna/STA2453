# 01_data_loading.R
# Retrieve file paths and load spatial transcriptomics data.
#
# This script extracts file paths for gene expression and metadata files,
# extracts sample names, and loads each dataset into list objects.
# It assumes that 'dir_path' is defined in 00_setup.R.
#
# Get file paths for expression and metadata files
expression_files <- dir_path[str_detect(dir_path, "expression")]
expression_files <- str_c(expression_files, "/", list.files(expression_files))

metadata_files <- dir_path[str_detect(dir_path, "metadata")]
metadata_files <- str_c(metadata_files, "/", list.files(metadata_files))

# Extract sample names from file names
sample_names <- str_remove(expression_files, fixed("~/Documents/***COURSEWORK***/STA2453H/Research/(+) GitHub/data/expression/datasets_mouse_brain_map_BrainReceptorShowcase_"))
sample_names <- str_remove(sample_names, "_cell_by_gene")
sample_names <- str_remove(sample_names, ".csv")

# Initialize list objects to store datasets (one element per sample)
data_list <- vector(mode = "list", length = length(expression_files))
names(data_list) <- sample_names
metadata_list <- vector(mode = "list", length = length(expression_files))
names(metadata_list) <- sample_names

# Load datasets: gene expression matrices and corresponding metadata
for (i in seq_along(expression_files)) {
  # Load gene expression matrix
  temp_expr <- as.data.frame(fread(expression_files[i]))[,-1]
  
  # Remove 'Blank' control probes from the gene list
  data_list[[i]] <- temp_expr[, str_detect(colnames(temp_expr), "Blank", negate = TRUE)]
  
  # Load metadata as a data frame
  temp_meta <- as.data.frame(fread(metadata_files[i]))
  colnames(temp_meta)[4:5] <- c("x","y")
  metadata_list[[i]] <- temp_meta
}
