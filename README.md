# Proteomics QC and Batch Correction Pipeline

## Overview

This repository contains an R-based pipeline for the preprocessing of quantitative proteomics data generated with Spectronaut.

The workflow is designed to perform:
- Data cleaning and filtering
- Quality control (PCA-based)
- Missing value assessment and imputation
- Batch effect correction using ComBat

This pipeline focuses on **data preprocessing** and **does not include differential expression analysis**. It generates a high-quality dataset ready for downstream statistical analysis.

---

## Workflow

### 1. Data Import
- Spectronaut output table is loaded
- Gene annotations are cleaned (removal of malformed entries)

### 2. Data Filtering
- Removal of:
  - Proteins without gene names
  - Low-confidence annotations (e.g. *putative*, *cDNA*, keratins)

### 3. Sample Processing
- Extraction of intensity columns
- Standardization of sample names
- Matching with experimental metadata

### 4. Data Structuring
- Creation of a `SummarizedExperiment` object using the DEP framework

### 5. Quality Control (Pre-filtering)
- Principal Component Analysis (PCA)
- Assessment of:
  - Biological condition
  - Batch effects
  - Technical variables

### 6. Missing Value Handling
- Filtering proteins with excessive missing values
- Visualization of missing data patterns

### 7. Imputation
- Missing values are imputed using K-Nearest Neighbors (KNN)

### 8. Batch Effect Correction
- Batch effects are corrected using **ComBat** (empirical Bayes method)
- Biological condition is preserved during correction

### 9. Post-processing QC
- PCA performed again to validate correction
- Detection of potential outliers

---


