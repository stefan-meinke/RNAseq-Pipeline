# ---------------------------------- #
# Script for generating PSI plots 
#                                                     
# Author: Dr. Stefan Meinke
# R version 4.3.0
#
# ---------------------------------- #

#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
  library(readr)
  library(readxl)
  library(magrittr)
  library(purrr)
  library(tidyr)
  library(GenomicRanges)
  library(rtracklayer)
  library(ggplot2)
  library(yaml)
  library(optparse)
})


# Define command-line options
option_list <- list(
  make_option(c("-d", "--psi_dir"), type="character", default=NULL,
              help="Directory containing psi files", metavar="dir"),
  make_option(c("-c", "--colData"), type="character", default=NULL,
              help="Path to the colData Excel file (sample metadata). Must include a 'Sample' column and at least one metadata column.", metavar="file"),
  make_option(c("-g", "--gtf"), type="character", default=NULL,
              help="Path to the GTF annotation file", metavar="file"),
  make_option(c("-t", "--target_genes"), type="character", default=NULL,
              help="Path to the target_genes.txt file containing one gene per row (no header)", metavar="file"),
  make_option(c("-o", "--output"), type="character", default="Results/PSI_plots",
              help="Output directory name for PSI plots [default %default]", metavar="file")
)

opt_parser <- OptionParser(option_list=option_list)
opts <- parse_args(opt_parser)

# Check for required arguments
if (is.null(opts$psi_dir)) { #|| is.null(opts$colData) || is.null(opts$gtf)
  print_help(opt_parser)
  stop("Error: --psi_dir, --colData, --gtf, and --target_genes are required.", call.=FALSE)
}

# Check if files/directories exist
if (!dir.exists(opts$psi_dir)) stop("Error: psi directory does not exist.")
if (!file.exists(opts$colData)) stop("Error: colData file does not exist.")
if (!file.exists(opts$gtf)) stop("Error: GTF file does not exist.")
if (!file.exists(opts$target_genes)) stop("Error: target_genes.txt fiel does not exist.")



# -------------------- #
# 1. Load the GTF file #
# -------------------- #
gtf_file <- opts$gtf
gtf <- rtracklayer::import(gtf_file)


# ----------------------------------------------- #
# 2. Load the sample information and target genes #
# ----------------------------------------------- #
colData <- read_xlsx(opts$colData)

config <- read_yaml("config.yml")

# Get the contrast definition
contrasts_list <- config$contrasts

target_genes <- readLines(opts$target_genes)

# ------------------------------------------------------- #
# 3. Load the PSI files and prepare the data for plotting #
# ------------------------------------------------------- #
psi_files <- opts$psi_dir
files_list_psi <- list.files(path = psi_files, pattern = "*.psi", full.names = TRUE, recursive = FALSE)
psi_dfs <- lapply(files_list_psi, function(x) read_tsv(x, col_names = c("exonID", "PSI", "low_inclusion_filter")))

names(psi_dfs) <- gsub(".psi","", basename(files_list_psi))

combined_list_psi <- list()
combined_list_psi <- imap(psi_dfs, ~ mutate(.x, Sample = .y))

PSI_all <- bind_rows(combined_list_psi) %>% 
  left_join(colData, by = "Sample")


# extract plot groups from config.yml
plot_groups_config <- config$plot_groups

# Get the full list of groups to plot and the reference group
plot_groups <- plot_groups_config$groups  
ref_level <- plot_groups_config$reference


# adjust PSI_all for plotting
PSI_all_adj <- PSI_all %>% 
  separate(exonID, into = c("geneID", "exon"), sep = ":") %>%
  mutate(exon = exon %>% as.numeric) %>% 
  mutate(group = factor(group, levels = c(ref_level, setdiff(plot_groups, ref_level))))




# -- 4. Load the psi_plot function --
psi_plot <- function(gene_name, psi_table, gtf, 
                     y_max = 115,
                     base_size = 10) {
  
  # 1) Filter GTF for exons of the given gene_name
  gene_exons <- gtf[gtf$type == "exon" & !is.na(gtf$gene_name) & gtf$gene_name == gene_name]
  if (length(gene_exons) == 0) {
    stop("No exons found for gene_name: ", gene_name)
  }
  
  # Convert GRanges to data frame and sort by genomic start
  exon_df <- as.data.frame(gene_exons) %>%
    arrange(start)
  
  # 2) Determine the strand from the gene-level record
  gene_gtf <- gtf[gtf$type == "gene" & !is.na(gtf$gene_name) & gtf$gene_name == gene_name]  
  if (length(gene_gtf) == 0) {
    stop("No gene-level record found for gene_name: ", gene_name)
  }
  
  strand_goi <- as.character(strand(gene_gtf))[1]  # assume single record
  
  # 3) filter psi table by gene of interest: extract corresponding gene_id
  goi_id <- unique(gene_gtf$gene_id)
  psi_table <- psi_table %>% filter(geneID == goi_id)
  
  # 4) Reorder exons if negative strand so that exon 1 is the 5′ end (leftmost in transcription order)
  if (strand_goi == "-") {
    psi_table$exon <- rev(psi_table$exon)
    exon_df <- exon_df %>% arrange(desc(start))
  } else {
    exon_df <- exon_df %>% arrange(start)
  }
  
  # Assign a new exon_index based on the sorted order
  exon_df <- exon_df %>% mutate(exon_index = row_number())
  
  # 4) Compute min/max genomic coords
  minCoord <- min(exon_df$start)
  maxCoord <- max(exon_df$end)
  genomic_length <- maxCoord - minCoord
  N <- nrow(exon_df)  # number of exons
  
  # 5) Define a scaling function so exons are drawn to scale at the top.
  #    For "+" strand, map minCoord -> 1, maxCoord -> N.
  #    For "-" strand, flip it so maxCoord -> 1, minCoord -> N.
  if (strand_goi == "+") {
    scale_genomic <- function(x) {
      1 + (N - 1) * (x - minCoord) / genomic_length
    }
  } else {
    scale_genomic <- function(x) {
      1 + (N - 1) * (maxCoord - x) / genomic_length
    }
  }
  
  # 6) Add scaled coordinates for each exon, then merge with your PSI table.
  #    Assume psi_table has columns: exon, PSI, Sample, Condition, etc.
  exon_df <- exon_df %>%
    mutate(
      scaled_start = scale_genomic(start),
      scaled_end   = scale_genomic(end),
      scaled_mid   = (scaled_start + scaled_end) / 2
    ) %>%
    left_join(psi_table, by = c("exon_index" = "exon")) %>%
    mutate(
      PSI    = PSI * 100,
      group = factor(group, levels = c(ref_level, setdiff(plot_groups, ref_level)))
    )
  
  # 7) Compute mean PSI for each exon & Condition (optional)
  meanPSI_df <- exon_df %>%
    group_by(exon_index, group) %>%
    summarize(
      mean_PSI = mean(PSI, na.rm = TRUE),
      sd_PSI   = sd(PSI, na.rm = TRUE),
      .groups  = "drop"
    )
  
  # 8) Attach the mean PSI back to the exon_df
  exon_df <- exon_df %>%
    left_join(meanPSI_df, by = c("exon_index", "group"))
  
  # 9) Build intron segments for the top track (scaled)
  intron_df <- exon_df %>%
    distinct(exon_index, scaled_start, scaled_end) %>%
    arrange(exon_index) %>%
    transmute(
      exon_index,
      intron_x1 = scaled_end,
      intron_x2 = lead(scaled_start)
    ) # %>%
    # slice(1:(nrow(.) - 1))
  
  # 10) Create partial vertical lines for the bottom (exon axis)
  vertical_lines_df <- data.frame(x = unique(exon_df$exon_index))
  
  # 11) Determine x-axis breaks. Use a sequence for breaks,
  #     and set x-axis limits to be c(0.5, N + 0.5)
  breaks_seq <- if (N > 50) seq(1, N, 5) else seq(1, N, 1)
  x_limits <- c(0.5, N + 0.5)
  
  # 12) Possibly rotate x-axis labels if many exons
  angle <- if (N > 50) 90 else 0
  hjust <- if (N > 50) 1 else 0.5
  
  # 13) Construct the ggplot
  p <- ggplot(exon_df) +
    theme_minimal(base_size = base_size) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks.y = element_line(color = "grey80"),
      axis.text.x = element_text(angle = angle, hjust = hjust),
      # axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      plot.margin = margin(2, 2, 2, 2, "pt")
    ) +
    
    # (A) Draw introns as horizontal lines at y = 112.5
    geom_segment(
      data = intron_df,
      aes(x = intron_x1, y = 112.5, xend = intron_x2, yend = 112.5),
      color = "grey10",
      linewidth = 0.1
    ) +
    
    # (B) Connector lines from y = 100 up to y = 110
    geom_segment(
      aes(x = exon_index, y = 100, xend = scaled_mid, yend = 110),
      color = "grey80",
      linewidth = 0.05
    ) +
    
    # (C) Draw the scaled exons as rectangles at the top
    geom_rect(
      aes(xmin = scaled_start, xmax = scaled_end,
          ymin = 110, ymax = 115),
      fill = "grey90", color = "grey10", linewidth = 0.1
    ) +
    
    # (D) Draw vertical lines from y = 0 to y = 100
    geom_segment(
      data = vertical_lines_df,
      aes(x = x, xend = x, y = 0, yend = 100),
      color = "grey80",
      linewidth = 0.1
    ) +
    
    # (E) Plot the mean PSI line for each Condition
    geom_line(aes(x = exon_index, 
                  y = mean_PSI, 
                  color = group, 
                  group = group)) +
    geom_ribbon(data = exon_df, aes(x = exon_index, ymin = mean_PSI - sd_PSI, ymax = mean_PSI + sd_PSI, fill = group), alpha = 0.2, color = NA) +
    scale_color_manual(values = c("grey60", "darkred", "#5B8FA8", "#ADB17D", "#FFB547")) +
    scale_fill_manual(values = c("grey60", "darkred", "#5B8FA8", "#ADB17D", "#FFB547")) +
    
    # (F) Scale and labels
    scale_x_continuous(
      breaks = breaks_seq,
      limits = x_limits,
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(0, y_max),
      breaks = c(0, 25, 50, 75, 100),
      expand = c(0, 0)
    ) +
    labs(
      x = "Exon",
      y = "PSI",
      title = bquote(italic(.(gene_name)))
    )
  
  return(p)
}


# Loop over each target gene, generate the plot, and save it as a PDF
for (gene in target_genes) {
  message("Plotting for gene: ", gene)
  
  # Compute the number of exons (N) for this gene from the GTF
  gene_exons <- gtf[gtf$type == "exon" & !is.na(gtf$gene_name) & gtf$gene_name == gene]
  if (length(gene_exons) == 0) {
    message("No exons found for gene: ", gene)
    next
  }
  # Convert to data frame and sort (also convert Rle columns if needed)
  gene_exon_df <- as.data.frame(gene_exons) %>%
    mutate(across(everything(), ~ if (inherits(.x, "Rle")) as.vector(.x) else .x)) %>%
    arrange(start)
  
  N <- nrow(gene_exon_df)
  
  # Calculate dynamic width using your formula:
  plot_width <- -2.47 + 2.3 * log(N)
  
  # (Optional: set a minimum or maximum width)
  plot_width <- max(3, min(plot_width, 11.5))
  
  # Generate the plot using your psi_plot function
  p <- psi_plot(gene, PSI_all_adj, gtf)
  
  # Define the output file name
  output_file <- file.path(opts$output, paste0("psiPlot_", gene, ".pdf")) 
  
  # Save the plot with the dynamically computed width and fixed height of 1.7 inches
  ggsave(filename = output_file, plot = p, width = plot_width, height = 1.7)
  
  message("Saved plot for gene: ", gene, " to: ", output_file)
}



