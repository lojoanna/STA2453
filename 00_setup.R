# 00_setup.R
# Load required libraries and set up project options.
#
# This script loads necessary libraries and defines global options,
# including specifying the base directory for the dataset.
#
library(tidyverse)
library(data.table)
library(Matrix)
library(Rfast)
library(matrixStats)
library(ggridges)
library(reticulate)
library(anndata)
library(ComplexHeatmap)
library(lmtest)

# Specify directory where data is saved
dir_path <- "~/Documents/***COURSEWORK***/STA2453H/Research/(+) GitHub/data"
dir_path <- str_c(dir_path, "/", list.files(dir_path))