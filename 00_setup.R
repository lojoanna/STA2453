# 00_setup.R
# Load required libraries and set up project options.
#
# This script loads necessary libraries and defines global options,
# including specifying the base directory for the dataset.
#

# Load required libraries

# Tidyverse + data handling
library(tidyverse)
library(data.table)
library(Matrix)
library(matrixStats)
library(Rfast)

# Plotting
library(ggridges)
library(ComplexHeatmap)
library(circlize)

# Modeling & statistics
library(lmtest)
library(mclust)
library(pROC)

# Parallelization
library(future.apply)

# Python interop
library(reticulate)
library(anndata)

# Additional utilities
library(conflicted)
library(hash)

# Specify directory where data is saved
dir_path <- "~/Documents/***COURSEWORK***/STA2453H/Research/(+) GitHub/data"
dir_path <- str_c(dir_path, "/", list.files(dir_path))