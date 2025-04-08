# ------------------------------------------------------------------- #
# Script for Differential Gene Expression (DGE) analysis using DESeq2 
#                                                     
# Author: Stefan Meinke
# R version 4.3.0
#
# ------------------------------------------------------------------- #

#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(tidyr)
  library(readxl)
  library(readr)
  library(dplyr)
  library(stringr)
  library(DESeq2)
  library(apeglm)
  library(purrr)
  library(GenomicFeatures)
  library(magrittr)
  library(writexl)
  library(tibble)
  library(yaml)
})

# Define command-line options
option_list <- list(
  make_option(c("-d", "--counts_dir"), type="character", default=NULL,
              help="Directory containing featureCounts output TXT files", metavar="dir"),
  make_option(c("-c", "--colData"), type="character", default=NULL,
              help="Path to the colData Excel file (sample metadata). Must include a 'Sample' column and at least one metadata column.", metavar="file"),
  make_option(c("-g", "--gtf"), type="character", default=NULL,
              help="Path to the GTF annotation file to extract geneID and gene name", metavar="file"),
  make_option(c("-o", "--output"), type="character", default="Results/DESeq2/DESeq_results.xlsx",
              help="Output file name for DESeq2 results [default %default]", metavar="file")
)

opt_parser <- OptionParser(option_list=option_list)
opts <- parse_args(opt_parser)

# Check for required arguments
if (is.null(opts$counts_dir) || is.null(opts$colData) || is.null(opts$gtf)) {
  print_help(opt_parser)
  stop("Error: --counts_dir, --colData, and --gtf are required.", call.=FALSE)
}

# Check if files/directories exist
if (!dir.exists(opts$counts_dir)) stop("Error: counts directory does not exist.")
if (!file.exists(opts$colData)) stop("Error: colData file does not exist.")
if (!file.exists(opts$gtf)) stop("Error: GTF file does not exist.")



# ----------------------------------------------------- #
# 1. Load and prepare the counts files for DGE analysis #
# ----------------------------------------------------- #

message("Loading featureCounts output files from: ", opts$counts_dir)

# List all files ending with "featureCounts.txt" (searching recursively)
reads <- list.files(
  path = opts$counts_dir,
  pattern = "*featureCounts.txt$",
  full.names = TRUE,
  recursive = TRUE
)
if (length(reads) == 0) {
  stop("Error: No featureCounts output files found in the specified directory.")
}

# Extract sample names from filenames by removing the suffix
samples <- gsub("_featureCounts.txt", "", basename(reads))

# Load and process each file in one step:
#  - Read the file (skipping the header line)
#  - Select the first column (gene ID) and the last column (counts)
#  - Rename the first column to "geneID" and the count column to the sample name
counts_list <- Map(function(file, sample) {
  df <- read_tsv(file, skip = 1, col_types = cols())
  first_col <- colnames(df)[1]
  last_col  <- colnames(df)[ncol(df)]
  df %>% 
    dplyr::select(geneID = !!rlang::sym(first_col), !!sample := !!rlang::sym(last_col))
}, reads, samples)


# Merge all dataframes by 'geneID'
counts <- reduce(counts_list, left_join, by = "geneID")


# Convert the merged dataframe to a counts matrix (with geneID as rownames)
counts_matrix <- counts %>% 
  as.data.frame() %>% 
  column_to_rownames("geneID")



# --------------------------------- #
# 2. Load colData (sample metadata) #
# --------------------------------- #

colData <- tryCatch({
  read_xlsx(opts$colData)
}, error = function(e) {
  stop("Error reading colData file: ", e$message)
})
if (!"Sample" %in% colnames(colData)) {
  stop("Error: colData file must contain a 'Sample' column.")
}

# Process colData metadata
metadata_cols <- setdiff(colnames(colData), "Sample")
if (length(metadata_cols) == 0) {
  stop("Error: colData file must have at least one metadata column besides 'Sample'.")
} else if (length(metadata_cols) == 1) {
  if (metadata_cols != "group") {
    message("Renaming metadata column '", metadata_cols, "' to 'group'.")
    colData <- colData %>% 
      dplyr::rename(group = !!metadata_cols)
  }
} else {
  # Combine multiple metadata columns into a single 'group' column (separated by underscores)
  colData <- colData %>%
    mutate(group = do.call(paste, c(across(
      all_of(metadata_cols)
    ), sep = "_")))
}

colData <- colData %>% column_to_rownames("Sample")


# --------------------------------------------- #
# 3. Extract gene annotations from the GTF file #
# --------------------------------------------- #

message("Extracting gene annotations from GTF file: ", opts$gtf)
txdb <- tryCatch({
  rtracklayer::import(opts$gtf)
}, error = function(e) {
  stop("Error loading GTF file: ", e$message)
})

gene_annotations <- as.data.frame(txdb) %>%
  dplyr::select(gene_id, gene_name) %>%
  dplyr::rename(geneID = gene_id, Gene = gene_name) %>% 
  drop_na %>% 
  unique
if(nrow(gene_annotations) == 0) {
  stop("Error: No gene annotations extracted from the GTF file.")
}



# ---------------------------------- #
# 5. Run DESeq2 and save the results #
# ---------------------------------- #

# make sure colnames(counts_matrix) and rownames(colData) are in the same order
common_samples <- intersect(colnames(counts_matrix), rownames(colData))
counts_matrix <- counts_matrix[, common_samples]
colData <- colData[common_samples, , drop = FALSE]


# Read config file to extract contrasts
config <- yaml::read_yaml("config.yml")
contrast_list <- config$contrasts

# Initialize a list to store results (optional)
results_list <- list()

# Loop over each contrast
for (contrast_name in names(contrast_list)) {
  contrast_vec <- contrast_list[[contrast_name]] 
  
  # Create a temporary copy of colData and relevel the 'group' factor so that the reference is the second element
  colData_temp <- colData
  colData_temp$group <- factor(colData_temp$group)  # ensure it's a factor
  colData_temp$group <- relevel(colData_temp$group, ref = contrast_vec[[2]])
  
  # Build DESeq2 dataset for this contrast
  dds_temp <- DESeqDataSetFromMatrix(countData = counts_matrix,
                                     colData = colData_temp,
                                     design = ~ group)
  
  dds_temp <- dds_temp[rowSums(counts(dds_temp)) >= 10, ]
  
  # Run DESeq2
  dds_temp <- DESeq(dds_temp)
  
  # Build the coefficient name.
  # With the factor releveling above, DESeq2 should create a coefficient of the form:
  # "group_<treatment>_vs_<reference>"
  coef_name <- paste0("group_", contrast_vec[[1]], "_vs_", contrast_vec[[2]])
  
  # Check that the coefficient is present
  rn <- resultsNames(dds_temp)
  if (!(coef_name %in% rn)) {
    stop("Coefficient ", coef_name, " not found. Check your factor levels and contrast definitions.")
  }
  
  # Perform LFC shrinkage using the coefficient name
  res_shrunk <- lfcShrink(dds_temp, coef = coef_name, type = "apeglm")
  res_df <- as.data.frame(res_shrunk) %>% rownames_to_column("geneID")
  
  # Join with gene annotations (assumes gene_annotations is defined)
  res_df <- res_df %>% 
    left_join(gene_annotations, by = "geneID") %>% 
    mutate(contrast = contrast_name)
  
  
  results_list[[contrast_name]] <- res_df
}


result_table <- bind_rows(results_list)

# Create output filename dynamically
output_file <- "Results/DESeq2/DESeq_results.xlsx"
  
# Write results to Excel
write_xlsx(result_table, output_file)
message("Differential expression analysis saved to: ", output_file)
  


