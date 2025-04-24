# ----------------------------------------- #
# Script for gene expression data analysis  
#                                                     
# Author: Stefan Meinke
# R version 4.3.0
#
# ----------------------------------------- #

#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(tidyr)
  library(readxl)
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(DESeq2)
  library(purrr)
  library(magrittr)
  library(writexl)
  library(tibble)
  library(yaml)
})

# Define command-line options
option_list <- list(
  make_option(c("-r", "--deseq_result"), type="character", default=NULL,
              help="Directory containing DESeq2 output excel file", metavar="dir"),
  make_option(c("-d", "--counts_dir"), type="character", default=NULL,
              help="Directory containing featureCounts output TXT files", metavar="dir"),
  make_option(c("-c", "--colData"), type="character", default=NULL,
              help="Path to the colData Excel file (sample metadata). Must include a 'Sample' column and at least one metadata column.", metavar="file"),
  make_option(c("-g", "--gtf"), type="character", default=NULL,
              help="Path to the GTF annotation file to extract geneID and gene name", metavar="file"),
  make_option(c("-f", "--config"), type="character", default=NULL,
              help="Path to the YAML config file specifying contrasts and filter cutoffs", metavar="file"),
  make_option(c("-o", "--output"), type="character", default="Results/DESeq2",
              help="Output folder name for enrichment results [default %default]", metavar="dir")
)

opt_parser <- OptionParser(option_list=option_list)
opts <- parse_args(opt_parser)

# Check for required arguments
if (is.null(opts$deseq_result) | is.null(opts$counts_dir) | is.null(opts$colData) | is.null(opts$gtf)) {
  print_help(opt_parser)
  stop("Error: --deseq_result, --counts_dir, --colData, and --gtf are required.", call.=FALSE)
}

# Check if files/directories exist
if (!dir.exists(opts$deseq_result)) stop("Error: DESeq2 result directory does not exist.")
if (!dir.exists(opts$counts_dir)) stop("Error: counts directory does not exist.")
if (!file.exists(opts$colData)) stop("Error: colData file does not exist.")
if (!file.exists(opts$gtf)) stop("Error: GTF file does not exist.")

# Define directory for outputs (relative to the Excel output)
output_dir <- file.path(dirname(opts$output))
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}



# ------------------- #
# custom ggplot theme #
# ------------------- #
my_theme <-   theme(line = element_line(color = "black"),
                    text = element_text(size = 10, color = "black"),
                    panel.background = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    strip.background = element_blank(),
                    axis.line = element_line(color = "black"),
                    strip.text.x = element_text(face = "bold"),
                    axis.text = element_text(color = "black", size = 10),
                    axis.ticks = element_line(color = "black"),
                    legend.key = element_rect(fill = NA, color = NA),
                    legend.position = "right",
                    legend.text = element_text(size = 10))

# --------------------------------- #
# Helper functions for saving plots #
# --------------------------------- #
# For ggplot2 objects
save_ggplot <- function(plot_object, file_path, width, height, dpi = 300) {
  out_dir <- dirname(file_path)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  ggsave(filename = file_path, plot = plot_object, width = width, height = height, dpi = dpi)
}


# ----------------------- #
# Load configuration file #
# ----------------------- #
message("Loading configuration from: ", opts$config)
config <- yaml::read_yaml(opts$config)


DGE_filter <- config$DGE_filter
padj_cutoff <- DGE_filter$padj
FC_cutoff <- DGE_filter$FC

# Convert the FC cutoff to log2 scale (for the log2FoldChange values)
log2_FC_cutoff <- log2(FC_cutoff)

contrast_list <- config$contrasts
contrast_vec <- contrast_list[[1]] 

# --------------------------- #
# Load the DESeq2 result file #
# --------------------------- #
message("Loading DESeq2 result file from: ", opts$deseq_result)

deseq_file <- file.path(opts$deseq_result, "DESeq_results.xlsx")
deseq <- read_xlsx(deseq_file)

deseq_filtered <- deseq %>% 
  filter(padj < padj_cutoff & abs(log2FoldChange) >= FC_cutoff)


# -------------------------------------------------- #
# Load and prepare the counts files for DGE analysis #
# -------------------------------------------------- #
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


# ------------------------------ #
# Load colData (sample metadata) #
# ------------------------------ #

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

# colData <- read_xlsx("colData.xlsx") # delete
colData <- colData %>% column_to_rownames("Sample")

# set the group as factor and define the reference level
colData$group <- factor(colData$group)

if (!(contrast_vec[[2]] %in% levels(colData$group))) {
  stop("Reference group in contrast_vec[[2]] not found in colData$group levels.")
}

colData$group <- relevel(colData$group, ref = contrast_vec[[2]])

# ------------------------------------------ #
# Extract gene annotations from the GTF file #
# ------------------------------------------ #

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


# ------------------------------------------------------------------------------ #
# Create DESeqDataSet, Pre-filter, and apply vst transformation for PCA analysis #
# ------------------------------------------------------------------------------ #
# make sure rownames(colData) = colnames(counts_matrix)
colData <- colData[colnames(counts_matrix), , drop = FALSE]

dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                              colData = colData,
                              design = ~ group) 

# filter out low count genes
dds <- dds[rowSums(counts(dds)) > 10, ]

# Use DESeq2's VST to transform counts to a scale where variance is more stable.
# Try vst() first, fallback to varianceStabilizingTransformation() if vst fails
vsd <- tryCatch({
  vst(dds, blind = FALSE)
}, error = function(e) {
  message("vst() failed: ", e$message)
  message("Falling back to varianceStabilizingTransformation()...")
  varianceStabilizingTransformation(dds, blind = FALSE)
})


expr_mat <- assay(vsd)


# -------------------- #
# Perform PCA analysis #
# -------------------- #

### Method 1: using prcomp()
# Run PCA on the expression matrix.
pca_res <- prcomp(t(expr_mat))

# Extract the percent variance explained.
pve <- (pca_res$sdev^2) / sum(pca_res$sdev^2) * 100

# Plot the PCA results.
pca_df <- data.frame(PC1 = pca_res$x[, 1],
                     PC2 = pca_res$x[, 2],
                     group = colData(dds)$group)


p1 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3, alpha = 0.8) +
  xlab(paste0("PC1: ", round(pve[1], 1), "% variance")) +
  ylab(paste0("PC2: ", round(pve[2], 1), "% variance")) +
  ggtitle("PCA on VST Counts (prcomp())") +
  scale_color_manual(values = c("grey70", "grey40", "grey10", "darkred", "#FFB547", "#ADB17D", "#5B8FA8")) +
  my_theme +  
  coord_fixed() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right")


save_ggplot(p1,
            file.path(output_dir, "PCA_prcomp.pdf"),
            width = 3.6,
            height = 3.0)


### Method 2: using plotPCA()
pcaData <- plotPCA(vsd, intgroup = c("group"), returnData = TRUE, ntop = 500)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# relevel the group factor
pcaData$group <- relevel(pcaData$group, ref = contrast_vec[[2]])

p2 <- ggplot(pcaData, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3, alpha = 0.8) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA on VST Counts (plotPCA())") +
  scale_color_manual(values = c("grey70", "grey40", "grey10", "darkred", "#FFB547", "#ADB17D", "#5B8FA8")) +
  my_theme +  
  coord_fixed() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right")

save_ggplot(p2,
            file.path(output_dir, "PCA_plotPCA.pdf"),
            width = 3.6,
            height = 3.0)


message("PCA analysis done!")


