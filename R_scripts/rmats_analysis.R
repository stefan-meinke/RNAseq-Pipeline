# -------------------------------------------------- #
# Script for the analysis of rMATS-turbo output data 
#                                                     
# Author: Dr. Stefan Meinke
# R version 4.3.0
#
# -------------------------------------------------- #

#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
  library(readr)
  library(writexl)
  library(readxl)
  library(magrittr)
  # library(ggbeeswarm)
  # library(ggsci)
  library(UpSetR)
  library(ggridges)
  library(scales)
  library(pheatmap)
  library(ggrepel)
  library(ggforce)
  library(optparse)
})

# Define command-line options
option_list <- list(
  make_option(c("-d", "--rmats_dir"), type="character", default=NULL,
              help="Directory containing rMATS output TXT files (*.MATS.JC.txt)", metavar="dir"),
  make_option(c("-o", "--output"), type="character", default="Results/rmats/rmats_results.xlsx",
              help="Output file name for rMATS results [default %default]", metavar="file")
)

opt_parser <- OptionParser(option_list=option_list)
opts <- parse_args(opt_parser)

# Check for required arguments
if (is.null(opts$rmats_dir)) {
  print_help(opt_parser)
  stop("Error: --rmats_dir", call.=FALSE)
}

# Check if files/directories exist
if (!dir.exists(opts$rmats_dir)) stop("Error: rMATS directory does not exist.")

# Define a dedicated directory for plot outputs (relative to the Excel output)
plot_output_dir <- file.path(dirname(opts$output), "plots")
if (!dir.exists(plot_output_dir)) {
  dir.create(plot_output_dir, recursive = TRUE)
}


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

# For base R plots (e.g. UpSetR output)
save_base_plot <- function(plot_expr, file_path, width, height) {
  out_dir <- dirname(file_path)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  pdf(file_path, width = width, height = height)
  eval(plot_expr)
  dev.off()
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

# ------------------------------------ #
# define color scheme for splice types #
# ------------------------------------ #
colors <- c("SE" = "#8F3931", 
            "MXE" = "#D6D6CE",
            "RI" = "#FFB547",
            "A3SS" = "#ADB17D",
            "A5SS" = "#5B8FA8")


# ------------------------------------------------------------------------- #
# utility function: sum comma-separated rMATS IncLevel values for filtering #
# ------------------------------------------------------------------------- #
sum_comma_separated <- function(x) {
  # Split the string by commas
  values <- strsplit(as.character(x), ",")[[1]]
  
  # Convert the resulting character vector to numeric
  numeric_values <- as.numeric(values)
  
  # Sum the numeric values
  sum_value <- sum(numeric_values, na.rm = TRUE)
  
  return(sum_value)
}



# ----------------------------------------------------------- #
# Function to load rMATS result files with optional filtering #
# ----------------------------------------------------------- #
load_rmats_data <- function(rmats_dir, filter_data = TRUE, fdr_threshold = 0.05, deltaPSI_threshold = 0.1) {
  # Initialize an empty list to store results
  rmats_results <- list()
  
  # List subdirectories within the base_dir
  subdirs <- list.dirs(rmats_dir, full.names = TRUE, recursive = FALSE)
  
  if (length(subdirs) == 0) {
    # If no subdirectories, process the base_dir itself
    subdirs <- rmats_dir
  }
  
  # Process each subdirectory
  for (subdir in subdirs) {
    # List all files in the subdirectory
    files_list <- list.files(
      path = subdir, 
      pattern = "*.MATS.JC.txt", 
      full.names = TRUE, 
      recursive = TRUE
    )
    
    # Read all files into a list of data frames
    dfs <- lapply(files_list, function(x) read_tsv(x, col_names = TRUE))
    
    # Extract the sample names and assign them as names to the list
    names(dfs) <- gsub("\\.MATS\\.JC\\.txt$", "", basename(files_list))
    
    # Optionally filter the data frames
    if (filter_data) {
      dfs_filtered <- lapply(dfs, function(df) {
        df %>%
          filter(
            sapply(IJC_SAMPLE_1, sum_comma_separated) >= 10 | 
              sapply(IJC_SAMPLE_2, sum_comma_separated) >= 10
          ) %>%
          filter(FDR < fdr_threshold & abs(IncLevelDifference) > deltaPSI_threshold) %>% 
          mutate(across(c(IJC_SAMPLE_1, IJC_SAMPLE_2, SJC_SAMPLE_1, SJC_SAMPLE_2), as.character))
      })
    } else {
      # If filter_data is FALSE, leave the data untouched.
      dfs_filtered <- dfs
    }
    
    # Add the sample name as a column and combine all data frames
    combined_list <- imap(dfs_filtered, ~ mutate(.x, Type = .y))
    AS_combined <- bind_rows(combined_list)
    
    if (nrow(AS_combined) == 0) {
      next
    }
    
    # Perform final processing: remove unneeded columns and rename key columns
    AS_combined %<>%
      dplyr::select(-c(ID...12, ID...14)) %>%
      dplyr::rename(ID = ID...1, geneID = GeneID, Gene = geneSymbol)
    
    # Use subdirectory name for naming list element
    subdir_name <- basename(subdir)
    rmats_results[[subdir_name]] <- AS_combined
  }
  
  return(rmats_results)
}



# ------------------------------------------- #
# 1. Load filtered and unfiltered rMATS files #
# ------------------------------------------- #
message("Loading rMATS-turbo output files from: ", opts$rmats_dir)

AS_all_list <- load_rmats_data(opts$rmats_dir, 
                               filter_data = TRUE, 
                               fdr_threshold = 0.01, 
                               deltaPSI_threshold = 0.15) 

# Combine the result files from each contrast into one dataframe
combined_list <- imap(AS_all_list, ~ mutate(.x, contrast = .y))



# Remove any data frames that are empty
non_empty_list <- Filter(function(df) nrow(df) > 0, combined_list)

if (length(non_empty_list) == 0) {
  message("All filtered data frames are empty. Nothing to process.")
  AS_all <- tibble()  # or exit/return if needed
} else {
  # Bind the non-empty data frames together
  AS_all <- bind_rows(non_empty_list)
  
  # Check for missing splice types
  expected_types <- c("SE", "MXE", "RI", "A3SS", "A5SS")
  present_types <- unique(AS_all$Type)
  missing_types <- setdiff(expected_types, present_types)
  if (length(missing_types) > 0) {
    message("The following splice types are missing after filtering: ", 
            paste(missing_types, collapse = ", "))
  }
  
  # Convert columns to factors as desired
  AS_all <- AS_all %>%  
    mutate(
      Type = factor(Type, levels = expected_types),
      contrast = factor(contrast)
    )
}


# save filtered result table
write_xlsx(AS_all, opts$output)


# load unfiltered splicing data
AS_all_unfiltered_list <- load_rmats_data(opts$rmats_dir, filter_data = FALSE) 
combined_list_unfiltered <- imap(AS_all_unfiltered_list, ~ mutate(.x, contrast = .y))
AS_all_unfiltered <- bind_rows(combined_list_unfiltered) %>% 
  mutate(Type = factor(Type, levels = c("SE", "MXE", "RI", "A3SS", "A5SS")),
         contrast = factor(contrast)
         )

  

# --------------------------------------------- #
# proportions of splice types for each contrast #
# --------------------------------------------- #

splicing_proportions <- AS_all %>% 
  group_by(contrast, Type) %>%
  dplyr::summarize(
    count = n(), # Count occurrences of each Type
    .groups = "drop" # Ungroup after summarizing
  ) %>%
  group_by(contrast) %>%
  mutate(
    proportion = count / sum(count) # Calculate proportion within each group
  ) %>%
  ungroup() %>% 
  mutate(Type = factor(Type, levels = c("A5SS", "A3SS", "RI", "MXE", "SE")))



p1 <- ggplot(splicing_proportions, aes(proportion, contrast, fill = Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors) +
  labs(y = "",
       fill = "") +
  my_theme +
  guides(fill = guide_legend(reverse = TRUE))


# dynamic height based on number of contrasts
n_contrasts <- length(unique(splicing_proportions$contrast))
base_height <- 1.5
height_p1 <- if (n_contrasts == 1) base_height else base_height + n_contrasts * 0.8

output_file_p1 <- file.path(plot_output_dir, "splice_proportions.pdf")
save_ggplot(p1, output_file_p1, width = 4, height = height_p1)
message("Saved splice proportions plot to: ", output_file_p1)


# ----------------------------------------------------------- #
# Number of splicing events for each contrast and Splice Type #
# ----------------------------------------------------------- #

n_events <- AS_all %>% 
  group_by(contrast, Type) %>% 
  dplyr::summarize(n = n()) %>% 
  mutate(Type = factor(Type, levels = c("A5SS", "A3SS", "RI", "MXE", "SE")))


p2 <- ggplot(n_events, aes(n, contrast, fill = Type)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = colors) +
  labs(x = "number of events",
       y = "") +
  facet_wrap(~ contrast, nrow = 1, scales = "free_x") +
  my_theme +
  guides(fill = guide_legend(reverse = TRUE))

base_height <- 1.8
height_p2 <- if (n_contrasts == 1) base_height else base_height + n_contrasts * 0.8
output_file_p2 <- file.path(plot_output_dir, "number_of_events.pdf")
save_ggplot(p2, output_file_p2, width = 4, height = height_p2)
message("Saved number of events plot to: ", output_file_p2)


# -------------------------------------------------------- #
# proportions of inclusion/exclusion events for each SF-KO #
# -------------------------------------------------------- #

inclusion_proportions <- AS_all %>% 
  mutate(splice_direction = ifelse(IncLevelDifference > 0, "inclusion", "exclusion")) %>% 
  group_by(contrast, splice_direction) %>% 
  dplyr::summarize(
    count = n(), # Count events in each splice direction
    .groups = "drop"
  ) %>%
  group_by(contrast) %>% 
  mutate(
    proportion = count / sum(count) # Calculate proportion within each group
  ) %>%
  ungroup() %>% 
  mutate(splice_direction = factor(splice_direction, levels = c("inclusion", "exclusion")))


p3 <- ggplot(inclusion_proportions, aes(proportion, contrast, fill = splice_direction)) +
  geom_bar(stat = "identity") +
  # scale_fill_uchicago("light") +
  scale_fill_manual(values = c("#D6604D", "#4393C3")) +
  labs(y = "",
       fill = "") +
  my_theme


base_height = 1
height_p3 <- if (n_contrasts == 1) base_height else base_height + n_contrasts * 0.8
output_file_p3 <- file.path(plot_output_dir, "inclusion_proportions.pdf")
save_ggplot(p3, output_file_p3, width = 4, height = height_p3)
message("Saved number of events plot to: ", output_file_p3)


# ----------------------------------------------- #
# overall deltaPSI for each splice Type and group #
# ----------------------------------------------- #

# label outliers
outliers <- AS_all %>%
  group_by(contrast, Type) %>%
  dplyr::summarize(
    Q1 = quantile(IncLevelDifference, 0.25),
    Q3 = quantile(IncLevelDifference, 0.75)
  ) %>%
  mutate(
    IQR = Q3 - Q1,
    LowerBound = Q1 - 1.5 * IQR,
    UpperBound = Q3 + 1.5 * IQR
  ) %>%
  left_join(AS_all, by = c("contrast", "Type")) %>%
  filter(IncLevelDifference < LowerBound | IncLevelDifference > UpperBound) %>%
  arrange(contrast, IncLevelDifference)


top_outliers <- outliers %>%
  group_by(contrast, Type) %>%
  top_n(n = 1, wt = IncLevelDifference) %>%
  bind_rows(
    outliers %>%
      group_by(contrast) %>%
      top_n(n = -1, wt = IncLevelDifference)
  ) %>%
  ungroup() %>% 
  dplyr::select(contrast, Q1, Q3, IQR, LowerBound, UpperBound, Gene, IncLevelDifference, Type)



p4 <- ggplot(AS_all, aes(IncLevelDifference, contrast, fill = Type)) +
  geom_boxplot() +
  labs(y = "",
       x = "delta PSI",
       fill = "") +
  geom_text(data = top_outliers, aes(label = Gene, x = IncLevelDifference), nudge_y = 0.05, check_overlap = TRUE, hjust = 0.5, vjust = -1, size = 3) +
  facet_wrap(~ contrast, nrow = 1) +
  # scale_fill_uchicago("light") +
  scale_fill_manual(values = colors) +
  scale_x_continuous(limits = c(-1.2, 1.2)) +
  my_theme




base_height <- 1.8
height_p4 <- if (n_contrasts == 1) base_height else base_height + n_contrasts * 0.8
output_file_p4 <- file.path(plot_output_dir, "deltaPSI_per_type.pdf")
save_ggplot(p4, output_file_p4, width = 4, height = height_p4)
message("Saved deltaPSI plot to: ", output_file_p4)




# --------------------------------------------------------------------------------- #
# Upset plot of overlapping genes among the splice types for all contrasts combined #
# --------------------------------------------------------------------------------- #
contrasts_vec <- as.character(unique(AS_all$contrast))
for (ctr in contrasts_vec) {
  tmp <- AS_all %>% 
    filter(contrast == ctr)
  
  upset_list <- list()
  for (stype in unique(tmp$Type)) {
    spliced_genes <- tmp %>% 
      filter(Type == stype) %>% 
      pull(Gene) %>% 
      unique()
    
    upset_list[[stype]] <- spliced_genes
  }
  
  output_file <- file.path(plot_output_dir, paste0("upsetPlot_", ctr, "_spliceTypes.pdf"))
  
  # Wrap the upset plot in a print() to ensure proper evaluation in a loop
  save_base_plot(expression(
    print(
      upset(fromList(upset_list),
            order.by = "freq",
            sets.x.label = "number of genes",
            text.scale = 1.5,
            point.size = 2,
            nintersects = NA,
            sets.bar.color = colors)
    )
  ), output_file, width = 5.6, height = 5.8)
  
  message("Saved upset plot for contrast ", ctr, " to: ", output_file)
}




# ----------------------------------------- #
# ranked dPSI plot of top 100 spliced genes #
# ----------------------------------------- #
spliceTypes_vec <- as.character(unique(AS_all$Type))

for (ctr in contrasts_vec) {
  tmp <- AS_all %>% filter(contrast == ctr)
  for (stype in spliceTypes_vec) {
    top50_inclusion <- tmp %>% 
      filter(Type == stype) %>% 
      arrange(desc(IncLevelDifference)) %>% 
      select(Gene, IncLevelDifference) %>% 
      drop_na() %>% 
      head(50)
    
    top50_exclusion <- tmp %>% 
      filter(Type == stype) %>% 
      arrange(IncLevelDifference) %>% 
      select(Gene, IncLevelDifference) %>% 
      drop_na() %>% 
      head(50)
    
    ranked_df <- bind_rows(top50_inclusion, top50_exclusion) %>% 
      arrange(desc(IncLevelDifference)) %>% 
      mutate(Gene = make.unique(Gene, sep = "."),
             Gene = factor(Gene, levels = Gene))
    
    p_ranked <- ggplot(ranked_df, aes(Gene, IncLevelDifference)) +
      geom_bar(stat = "identity", fill = "#ADB17D", color = "black") +
      labs(x = "", y = "dPSI [%]",
           title = paste0("Top ", stype, " events for contrast ", ctr)) +
      my_theme +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    
    output_file <- file.path(plot_output_dir, paste0("rankedPSI_", ctr, "_", stype,".pdf"))
    save_ggplot(p_ranked, output_file, width = 12, height = 2.2)
    
    message("Saved ranked PSI plot for contrast ", ctr, " and splice type ", stype, " to: ", output_file)
  }
}



# ------------------- #
# Generate PSI matrix #
# ------------------- #
# 
# for(type in spliceTypes){
#   
#   tmp <- AS_all %>% 
#     filter(Type == type) %>% 
#     dplyr::select(ID, Gene, IncLevel1, IncLevel2, contrast) %>% 
#     separate_rows(IncLevel1, IncLevel2, sep = ",")
# }
# 
# 
# tmp <- AS_all %>% 
#   filter(Type == "SE") %>% 
#   dplyr::select(ID, Gene, IncLevel1, IncLevel2, contrast) %>% 
#   separate_rows(IncLevel1, IncLevel2, sep = ",") %>% 
#   pivot_longer(cols = matches("IncLevel"), names_to = "group", values_to = "PSI") %>% 
#   group_by(contrast, ID, Gene, group) %>% 
#   mutate(Sample = paste0(group, "_", row_number())) %>% 
#   ungroup
#   
# 
# psi_mat <- tmp %>%
#   filter(!is.na(Gene)) %>% 
#   mutate(event = paste0(ID,":", Gene)) %>% 
#   dplyr::select(event, Sample, PSI) %>% 
#   pivot_wider(names_from = "Sample", values_from = "PSI") %>%
#   drop_na %>% 
#   column_to_rownames("event") %>% 
#   mutate(across(everything(), as.numeric))
# 
# 
# 
# 
# # -------------------------------- #
# # Heamtap of top 100 spliced genes #
# # -------------------------------- #
# 
# heatmap_list <- list()
# for(type in spliceTypes){
#   
#   top50_inclusion <- AS_all %>% 
#     filter(Type == type) %>% 
#     arrange(desc(IncLevelDifference)) %>% 
#     filter(!is.na(Gene)) %>%
#     mutate(event = paste0(ID,":",Gene)) %>% 
#     dplyr::select(event, IncLevelDifference) %>% 
#     head(50) 
#   
#   top50_exclusion <- AS_all %>% 
#     filter(Type == type) %>% 
#     arrange(IncLevelDifference) %>% 
#     filter(!is.na(Gene)) %>%
#     mutate(event = paste0(ID,":",Gene)) %>% 
#     dplyr::select(event, IncLevelDifference) %>% 
#     head(50) 
#   
#   top100_df <- bind_rows(top50_inclusion, top50_exclusion) 
#   
#   psi_mat <- AS_all %>% 
#     filter(Type == type,
#            !is.na(Gene)) %>% 
#     mutate(event = paste0(ID,":", Gene)) %>% 
#     dplyr::select(event, IncLevel1, IncLevel2, contrast) %>% 
#     separate_rows(IncLevel1, IncLevel2, sep = ",") %>% 
#     pivot_longer(cols = matches("IncLevel"), names_to = "group", values_to = "PSI") %>% 
#     group_by(contrast, event, group) %>% 
#     mutate(Sample = paste0(group, "_", row_number())) %>% 
#     ungroup %>% 
#     dplyr::select(event, Sample, PSI) %>% 
#     pivot_wider(names_from = "Sample", values_from = "PSI") %>%
#     drop_na %>% 
#     column_to_rownames("event") %>% 
#     mutate(across(everything(), as.numeric))
#   
#   psi_mat_top100 <- psi_mat[rownames(psi_mat) %in% top100_df$event, ]
#   
#   inc1_cols <- grep("IncLevel1", colnames(psi_mat_top100), value = TRUE)
#   inc2_cols <- grep("IncLevel2", colnames(psi_mat_top100), value = TRUE)
#   ordered_cols <- c(inc1_cols, inc2_cols)
#   psi_mat_top100 <- psi_mat_top100[, ordered_cols]
#   
#   
#   p <- pheatmap(psi_mat_top100,
#            # scale = "row",
#            cluster_rows = TRUE,
#            cluster_cols = FALSE,
#            show_rownames = TRUE,
#            show_colnames = TRUE,
#            treeheight_row=0, 
#            treeheight_col=0,
#            color = colorRampPalette(c("#4393C3", "white", "#D6604D"))(20),
#            main = type)
# 
#   heatmap_list[[type]] <- p
# }
# 
# 
# 
# 
# # ------------ #
# # Volcano Plot #
# # ------------ #
# 
# volcano_list <- list()
# 
# for (c in contrasts) {
#   # Create an empty nested list for each contrast
#   volcano_list[[c]] <- list()
#   
#   for (splice_type in spliceTypes) {
#     # Filter the data frame for the current contrast and splice type
#     volcano_df <- AS_all %>% 
#       filter(contrast == c,     # adjust if your contrast column has a different name
#              Type == splice_type, 
#              !is.na(Gene)) %>% 
#       mutate(event = paste0(ID, ":", Gene)) %>% 
#       select(event, IncLevelDifference, FDR) %>% 
#       mutate(FDR = ifelse(FDR == 0, min(FDR[FDR > 0], na.rm = TRUE), FDR),
#              splicing = ifelse(IncLevelDifference > 0, "inclusion", "exclusion"))
#     
#     # Generate the volcano plot
#     p <- ggplot(volcano_df, aes(IncLevelDifference, -log10(FDR), color = splicing)) +
#           geom_point() +
#           geom_vline(xintercept = c(-0.15, 0.15), linetype ="dashed") + # set lines at deltaPSI limit (based on your filter) 
#           geom_hline(yintercept = -log10(0.01), linetype = "dashed") + # set lines at FDR filter
#           geom_text_repel(data = volcano_df %>% filter(splicing =="inclusion") %>% slice_min(FDR, n = 15), max.overlaps = 15, aes(label = `event`), size = 2.5) +
#           geom_text_repel(data = volcano_df %>% filter(splicing == "exclusion") %>% slice_min(FDR, n = 15), max.overlaps = 15, aes(label = `event`), size = 2.5) +
#           scale_color_manual(values = c("inclusion" = "#4393C3", "exclusion" = "#D6604D")) +
#           labs(x = "deltaPSI",
#                y = "-log10(FDR)") +
#           # scale_x_continuous(limits=c(-6,6,1)) +
#           # scale_y_continuous(limits=c(NA, 600)) +
#           my_theme +
#           theme(legend.position = "none")
#     
#     # Save the plot in the list, using contrast and splice_type as keys
#     volcano_list[[c]][[splice_type]] <- p
#   }
# }
# 
# 
# volcano_list$KO_vs_WT$A5SS
# 
# 
# 
# 
# 
# # ------------------- #
# # Enrichment analysis #
# # ------------------- #
# 
# library(clusterProfiler)
# library(org.Hs.eg.db) # change to org.Mm.eg.db for mouse data
# 
# ################################
# ### KEGG enrichment analysis ###
# ################################
# 
# # analysis for each contrast across all splice types (taking all significant splicing events into account)
# background_genes <- unique(AS_all_unfiltered$Gene)
# 
# background_ids <- bitr(background_genes,
#                        fromType = "SYMBOL",
#                        toType = "ENTREZID",
#                        OrgDb = org.Hs.eg.db) #change to org.Mm.eg.db for mouse data
# KEGG_plots <- list()
# KEGG_dfs <- list()
# for(c in contrasts){
#   sign_genes <- AS_all %>% 
#     filter(contrast == c) %>% 
#     pull(Gene) %>% 
#     unique
#   
#   entrez_ids <- bitr(sign_genes, 
#                      fromType = "SYMBOL",
#                      toType = "ENTREZID",
#                      OrgDb = org.Hs.eg.db)
#   
#   kegg <- enrichKEGG(gene = entrez_ids$ENTREZID,
#                      organism = "hsa", # change to "mmu" for mouse data
#                      pvalueCutoff = 0.05,
#                      pAdjustMethod = "BH",
#                      universe = background_ids$ENTREZID)
#   
#   kegg_df <- as.data.frame(kegg) %>%
#     # mutate(Description = gsub(" - Mus musculus \\(house mouse\\)", "", Description)) %>% # include for mouse data
#     arrange(p.adjust) %>% 
#     mutate(Description = factor(Description, levels = rev(unique(Description))),
#            padj = ifelse(`p.adjust` < 0.05, "padj < 0.05", "padj > 0.05"),
#            `-log10(padj)` = -log10(p.adjust))
#   
#   p <- ggplot(head(kegg_df, 25), aes(x = `-log10(padj)`, y = Description)) +
#     geom_bar(stat = "identity", fill = "grey90", color = "black") +
#     geom_text(aes(label = Count), vjust = 0.5, hjust = -0.5) +
#     geom_vline(xintercept = -log10(0.05), linetype = "dashed") +
#     scale_x_continuous(limits = c(0, max(kegg_df$`-log10(padj)` + .5))) +
#     labs(x = "-log10(padj)",
#          y = "Enriched Term",
#          size = "#Genes",
#          fill = "") +
#     ggtitle(paste0("KEGG: ", c)) +
#     my_theme + 
#     theme(legend.title = element_text(),
#           legend.position = c(0.7, 0.15),
#           panel.border = element_blank(),
#           panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank())
#   
#   KEGG_plots[[c]] <- p
#   KEGG_dfs[[c]] <- kegg_df
# }
# 
# 
# KEGG_plots$KO_vs_WT
# 
# KEGG_dfs$KO_vs_WT
#   
# 
# #save the plot as pdf file
# # for(p in seq_along(KEGG_plots)){
# #   plot_name <- paste0(dir, paste0(format(Sys.Date(), "%Y%m%d"), "_KEGG_enrichment_", names(KEGG_plots)[p], ".pdf") # adapt dir as needed
# #   ggsave(plot_name, plot = KEGG_plots[[p]], width = 8, height = 6) # adjust width and height based on your needs
# # }
# 
# 
# 
# 
# ######################################
# ### Perform GO Enrichment Analysis ###
# ######################################
# 
# sign_genes <- AS_all %>%
#   filter(Type == "SE", # adapt all filter as needed 
#          contrast == "KO_vs_WT") %>% 
#   pull(Gene)
# 
# entrez_ids <- bitr(sign_genes, 
#                    fromType = "SYMBOL",
#                    toType = "ENTREZID",
#                    OrgDb = org.Hs.eg.db)
# 
# background_genes <- unique(AS_all_unfiltered$Gene)
# 
# background_ids <- bitr(background_genes,
#                        fromType = "SYMBOL",
#                        toType = "ENTREZID",
#                        OrgDb = org.Hs.eg.db)
# 
# 
# # define the function to perform the GO enrichment, generate the result dataframe and make the plot
# Go_plot <- function(genes,background,ont){ # genes = entrez ids of significant genes, background = entrez ids of background genes, ont = "BP" or "MF" or "CC"
#   go <- enrichGO(gene = genes,
#                  universe = background,
#                  OrgDb = org.Hs.eg.db, # or org.Mm.eg.db for mouse
#                  ont = ont,
#                  pvalueCutoff = 0.05,
#                  pAdjustMethod = "BH") # adjustment using Benjamini-Hochberg method, for other methods see ?enrichGO
#   
#   go_df <- go@result %>% 
#     dplyr::select(Description, p.adjust, ENTREZID = geneID, Count) %>% 
#     arrange(p.adjust) %>% 
#     filter(p.adjust < 0.05) %>% 
#     separate_rows(ENTREZID, sep = "/") %>% 
#     left_join(entrez_ids, by = "ENTREZID") %>% 
#     dplyr::rename(Gene = SYMBOL) %>% 
#     mutate(ont = ont)
#   
#   enriched_terms <- go$Description
#   p_values <- go$pvalue
#   padj <- go$p.adjust
#   gene_counts <- go$Count
#   
#   result_df <- data.frame(
#     Enriched_Term = enriched_terms,
#     P_Value = p_values,
#     Adjusted_P_Value = padj,
#     Gene_Count = gene_counts)
#   
#   result_df %<>% 
#     arrange(P_Value)
#   
#   top20 <- head(result_df, 20) # adjust for more/less pathways to be shown
#   
#   # Reorder the Enriched_Term factor variable based on the P_Value column
#   top20$Enriched_Term <- factor(top20$Enriched_Term, levels = top20$Enriched_Term)
#   
#   # Create the horizontal bar plot
#   plot <- ggplot(top20, aes(x = Adjusted_P_Value, y = reorder(Enriched_Term, -Adjusted_P_Value))) + # x = P_Value
#     # geom_point(shape=21, aes(size=Gene_Count, fill = Adjusted_P_Value)) +
#     geom_bar(stat = "identity", fill = "grey90", color = "black") +
#     # geom_text(aes(label = Gene_Count), vjust = 0.5, hjust=-0.5) +
#     labs(x = "adjusted P-Value", y = "")+ #,
#     ggtitle(if(ont == "BP"){
#       "GO Biological process"
#     } else {
#       if(ont == "MF"){
#         "GO Molecular function"
#       } else {
#         if(ont == "CC"){
#           "GO Cellular Component"
#         }}})+
#     geom_vline(xintercept = 0.05, linetype = "dashed") +
#     scale_x_continuous(trans = trans_reverser("log10")) +
#     # scale_size_continuous(range = c(3, 8)) +
#     my_theme +
#     theme(legend.title = element_text(),
#           legend.position = c(1, 0), legend.justification = c(1, 0)) # adjust if necessary
#   
#   return(list(plot = plot, go_df = go_df))
#   
# }
# 
# 
# # example for creating the GO Biological Process enrichment
# GO_plots <- list()
# GO_dfs <- list()
# for (ont in c("BP", "MF", "CC")) {
#   result <- Go_plot(entrez_ids$ENTREZID, background_ids$ENTREZID, ont)
#   
#   # Store plot and dataframe in respective lists
#   GO_plots[[ont]] <- result$plot
#   GO_dfs[[ont]] <- result$go_df
# }
# 
# 
# GO_plots$CC
# GO_dfs$CC


#save the plot as pdf file
# for(p in seq_along(GO_plots)){
#   plot_name <- paste0(dir, paste0(format(Sys.Date(), "%Y%m%d"), "_GO_enrichment_", names(GO_plots)[p], ".pdf")
#   ggsave(plot_name, plot = GO_plots[[p]], width = 8.3, height = 5.5)
# }