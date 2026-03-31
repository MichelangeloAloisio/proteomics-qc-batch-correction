# ==============================================================
# PROTEOMICS DATA ANALYSIS PIPELINE (DEP + ComBat)
# ==============================================================



# ==============================================================
# 1. Clean environment
# ==============================================================
rm(list = ls())
options(scipen = 9999999)


# ==============================================================
# 2. Set working directory
# ==============================================================
#setwd("D:/Desktop/progetti_v1/34_PROTEOMICA_METABOLOMICA_IBS/risultati_blendi_R")
setwd('D:/Desktop/progetti_v1/34_PROTEOMICA_METABOLOMICA_IBS/risultati_BLENDI_ROC')


# ==============================================================
# 3. Load required libraries
# ==============================================================
library(DEP)
library(tidyverse)
library(SummarizedExperiment)
library(sva)



# ==============================================================
# SECTION A — CREATE SummarizedExperiment OBJECT
# ==============================================================


# ==============================================================
# 4. Load proteomics data (Spectronaut export)
# ==============================================================
#data_raw <- read_delim("Cartella_merge_Spectranout_SHORT.csv", delim = ";", show_col_types = FALSE)

# Load Spectronaut output table
data_raw <- read_delim(
  "00_files_input_spectronout/spectronout_Libreria_20.000.csv",
  delim = ",",
  show_col_types = FALSE
)

# Remove leading/trailing semicolons from gene column
data_raw <- data_raw %>%
  mutate(PG.Genes = gsub("^;+|;+$", "", PG.Genes, perl = TRUE))

dim(data_raw)  # Expected: 9730 x 42


# ==============================================================
# 5. Data filtering
# ==============================================================

## FILTER 1: Remove rows with missing gene annotation
data_raw <- data_raw[!is.na(data_raw$PG.Genes), ]
dim(data_raw)

## FILTER 2: Remove low-confidence or unwanted protein annotations
pattern <- "putative|probable|cDNA|KERATINE|SIMILAR"
data_raw <- data_raw[
  !grepl(pattern, data_raw$PG.ProteinDescriptions, ignore.case = TRUE),
]
dim(data_raw)


# ==============================================================
# 6. Select and rename intensity columns
# ==============================================================

# Keep identifiers + intensity columns
data_filtered <- data_raw %>%
  select(
    PG.ProteinAccessions,
    PG.Genes,
    matches("\\[\\d+\\] \\d+[A-Z]\\.raw\\.PG\\.Quantity")
  )

dim(data_filtered)

# Extract intensity column names
quantity_cols <- names(data_filtered)[
  !names(data_filtered) %in% c("PG.ProteinAccessions", "PG.Genes")
]

# Rename columns to standardized format (e.g., s_A1)
new_names <- sapply(quantity_cols, function(x) {
  match <- regmatches(x, regexpr("\\d+[A-Z]", x))
  if (length(match) == 0) return(x)
  
  num <- gsub("[A-Z]", "", match)
  letter <- gsub("\\d+", "", match)
  
  paste0("s_", letter, num)
})

names(data_filtered)[names(data_filtered) %in% quantity_cols] <- new_names

data_final <- data_filtered
dim(data_final)


# ==============================================================
# 7. Create unique protein identifiers
# ==============================================================

data_unique <- make_unique(
  data_final,
  "PG.Genes",
  "PG.ProteinAccessions",
  delim = ";"
)

dim(data_unique)

# Identify intensity columns
intensity_cols <- which(startsWith(colnames(data_unique), "s_A"))


# ==============================================================
# 8. Load experimental design (metadata)
# ==============================================================

expdesign <- read_delim(
  "00_files_input_spectronout/metadata.csv",
  delim = ";",
  show_col_types = FALSE
)

dim(expdesign)

# Check consistency between metadata and data
all(expdesign$label %in% colnames(data_final))


# ==============================================================
# 9. Check for zero values before log2 transformation
# ==============================================================

dataframe_pre_log2_tranformation <- data_unique[
  , !colnames(data_unique) %in% c("name", "ID", "PG.ProteinAccessions", "PG.Genes")
]

has_zeros <- any(dataframe_pre_log2_tranformation == 0, na.rm = TRUE)

if (has_zeros) {
  cat("Warning: zero values detected (will become -Inf after log2).\n")
} else {
  cat("No zero values detected. Safe for log2 transformation.\n")
}


# ==============================================================
# 10. Create SummarizedExperiment object
# ==============================================================

data_se <- make_se(
  proteins_unique = data_unique,
  columns = intensity_cols,
  expdesign = expdesign
)

dim(data_se)


# ==============================================================
# 11. Initial PCA diagnostics
# ==============================================================

plot_pca(data_se, indicate = "condition", n = 6000)
plot_pca(data_se, indicate = "batch", n = 6000)
plot_pca(data_se, indicate = "breadford", n = 6000)


sum(is.na(assay(data_se)))


# ==============================================================
# 12. Filter proteins based on missing values
# ==============================================================

plot_frequency(data_se)

data_filt <- filter_missval(data_se, thr = 0)

dim(data_filt)

plot_frequency(data_filt)
plot_numbers(data_filt)
plot_detect(data_filt)

sum(is.na(assay(data_filt)))

plot_pca(data_filt, indicate = "batch", n = 6000)
plot_pca(data_filt, indicate = "condition", n = 6000)


# ==============================================================
# 13. Missing value imputation (KNN)
# ==============================================================

plot_missval(data_filt)

data_imp <- impute(data_filt, fun = "knn", rowmax = 0.1)

plot_imputation(data_filt, data_imp)

dim(data_imp)
sum(is.na(assay(data_imp)))

# PCA after imputation
plot_pca(data_imp, indicate = "batch", n = 6000)
plot_pca(data_imp, indicate = "condition", n = 6000)


# ==============================================================
# 14. Batch effect correction (ComBat)
# ==============================================================

# Prepare phenotype data
pheno <- as.data.frame(colData(data_imp))

# Extract expression matrix (log2 scale)
edata <- assay(data_imp)

# Check alignment
all(rownames(pheno) %in% colnames(edata))

# Define batch and model
batch <- pheno$batch
mod <- model.matrix(~ as.factor(condition), data = pheno)

# Apply ComBat
combat_edata3 <- ComBat(
  dat = edata,
  batch = batch,
  mod = mod,
  par.prior = TRUE,
  ref.batch = 3
)


# =======================================================================================
# SECTION B — REBUILD SummarizedExperiment AFTER ComBat to perform PCA control analysis
# ======================================================================================

# Convert back to linear scale
combat_edata3_DF <- as.data.frame(combat_edata3)
combat_edata3_DF <- 2^(combat_edata3_DF)

# Add identifiers
combat_edata3_DF$name <- rownames(combat_edata3_DF)
combat_edata3_DF$ID <- rownames(combat_edata3_DF)

# Identify intensity columns
intensity_cols_combat_edata3_DF <- which(
  startsWith(colnames(combat_edata3_DF), "SANA") |
    startsWith(colnames(combat_edata3_DF), "ALTERATA")
)

# Rename columns using metadata labels
new_colnames <- expdesign$label
colnames(combat_edata3_DF)[1:30] <- new_colnames

# Ensure metadata consistency
expdesign$sample <- expdesign$label

# Recreate SummarizedExperiment
data_se <- make_se(
  proteins_unique = combat_edata3_DF,
  columns = intensity_cols_combat_edata3_DF,
  expdesign = expdesign
)


# ==============================================================
# 15. Final PCA after batch correction
# ==============================================================

plot_pca(data_se, indicate = "batch", n = 6000)
plot_pca(data_se, indicate = "condition", n = 6000)


# ==============================================================
# 16. Save combat batch corrected file
# ==============================================================

write_csv(combat_edata3_DF, "combat_corrected_counts.csv")
