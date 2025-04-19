# Predicting Gene Expression from Neighbouring Cells

## Overview

This repository contains code for an independent research project analyzing spatial transcriptomics data from the Vizgen MERSCOPE Mouse Brain Map. The study investigates spatial patterns of gene expression, cell-neighbor relationships, and co-expression networks across multiple mouse brain slices.

The project applies preprocessing, exploratory data analysis (EDA), and statistical modeling to understand how spatial organization influences gene expression. The analysis leverages high-resolution spatial transcriptomic data, where individual cells' gene expression profiles are mapped to their physical locations.

## Data Source

The raw datasets are publicly available from **Vizgen's Mouse Brain Map**. Due to file size limitations, the data files are not included in this repository, but they can be accessed at:

📌 [Vizgen Mouse Brain Map](https://info.vizgen.com/mouse-brain-map)

To reproduce the analysis, follow these steps:

1. Click the link above to visit the Vizgen Mouse Brain Map webpage.
2. Scroll down to **Section 3: Download the Data**.
3. For **each brain slice**, download the following two files:
   - `cell_by_gene_SxRy.csv` — the gene expression matrix (cells × genes)
   - `cell_metadata_SxRy.csv` — the metadata for each cell
4. Here, `Sx` refers to the **slice number** (`S1`, `S2`, `S3`) and `Ry` refers to the **replicate** (`R1`, `R2`, `R3`), giving a total of **9 unique slice-replicate combinations** (e.g., `S1R1`, `S1R2`, ..., `S3R3`).
5. This results in **18 CSV files** in total (9 slices × 2 files per slice).
6. Save these files in a local directory and update the relevant **file paths** in `01_data_loading.R` to match your system.

> ⚠️ Do not rename the files. The script relies on the naming pattern `cell_by_gene_SxRy.csv` and `cell_metadata_SxRy.csv` to correctly identify and load each slice.

## Repository Structure

The repository is structured into modular scripts to ensure reproducibility and organization:
```
├── 00_setup.R                 # Loads libraries and sets global options
├── 01_data_loading.R          # Retrieves file paths and loads datasets
├── 02_analysis_functions.R    # Defines custom functions for analysis
├── 03_preprocessing.R         # Performs data cleaning and transformations
├── 04_exploratory_analysis.R  # Conducts exploratory data analysis (EDA)
├── 05_statistical_analysis.R  # Statistical modeling & evaluation
└── README.md                  # Project documentation
```

### Description of Scripts

- **`00_setup.R`** — Loads required libraries and sets global options.
- **`01_data_loading.R`** — Extracts file paths and loads gene expression and metadata.
- **`02_analysis_functions.R`** — Defines reusable functions for coordinate transformation, visualization, and analysis.
- **`03_preprocessing.R`** — Cleans the data by removing low-quality cells and rotating/normalizing spatial coordinates.
- **`04_exploratory_analysis.R`** — Visualizes total RNA counts, cell volume, gene-gene co-expression, and spatial correlations.
- **`05_statistical_analysis.R`** — Contains statistical modeling, hypothesis testing, and evaluation metrics.

## Reproducibility

To run the analysis:
1. Clone this repository:
   ```bash
   git clone https://github.com/yourusername/spatial-transcriptomics.git
   cd spatial-transcriptomics
   ```
2. Install the required R packages:
   ```r
   install.packages(c("tidyverse", "data.table", "Matrix", "Rfast", "matrixStats", 
                   "ggridges", "reticulate", "anndata", "ComplexHeatmap", "lmtest"))
   ```
4. Download the dataset from Vizgen Mouse Brain Map and place the files in your local directory.
5. Run the scripts in order:
   ```r
   source("00_setup.R")
   source("01_data_loading.R")
   source("02_analysis_functions.R")
   source("03_preprocessing.R")
   source("04_exploratory_analysis.R")
   source("05_statistical_analysis.R")
   ```
