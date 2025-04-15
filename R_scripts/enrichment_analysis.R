# ------------------------------------------------------------------- #
# Script for enrichment analysis of differentially expressed genes  
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
  library(dplyr)
  library(stringr)
  library(clusterProfiler)
  library(magrittr)
  library(writexl)
  library(ggplot2)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(ggforce)
  library(tibble)
  library(yaml)
})


# Define command-line options
option_list <- list(
  make_option(c("-d", "--deseq_result"), type="character", default=NULL,
              help="Directory containing DESeq2 output excel file", metavar="dir"),
  make_option(c("-a", "--annotation"), type = "character", default=NULL,
              help="Specify the organism ('m' for mouse, 'h' for human).", metavar="file"),
  make_option(c("-c", "--config"), type = "character", default="config.yml",
              help="Path to config file [default %default]", metavar="file"),
  make_option(c("-o", "--output"), type="character", default="Results/DESeq2/Enrichment",
              help="Output folder name for enrichment results [default %default]", metavar="dir")
)

opt_parser <- OptionParser(option_list=option_list)
opts <- parse_args(opt_parser)

# Check for required arguments
if (is.null(opts$deseq_result) | is.null(opts$annotation)) {
  print_help(opt_parser)
  stop("Error: --deseq_result and --annotation are required.", call.=FALSE)
}

# Check if files/directories exist
if (!dir.exists(opts$deseq_result)) stop("Error: DESeq2 result directory does not exist.")

# Define directory for outputs (relative to the Excel output)
output_dir <- file.path(dirname(opts$output))
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

output_dir <- "Results/DESeq2/Enrichment"
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


# -------------------------- #
# 1. Load configuration file #
# -------------------------- #
message("Loading configuration from: ", opts$config)
config <- yaml::read_yaml(opts$config)

config <- yaml::read_yaml("config.yml") ### delete


DGE_filter <- config$DGE_filter
padj_cutoff <- DGE_filter$padj
FC_cutoff <- DGE_filter$FC

# Convert the FC cutoff to log2 scale (for the log2FoldChange values)
log2_FC_cutoff <- log2(FC_cutoff)

# ------------------------------------------ #
# 2. Load and prepare the DESeq2 result file #
# ------------------------------------------ #
message("Loading DESeq2 result file from: ", opts$deseq_result)
deseq_file <- file.path(opts$deseq_result, "DESeq_results.xlsx")
deseq <- read_xlsx(deseq_file)


deseq <- read_xlsx("Results/DESeq2/DESeq_results.xlsx") ### delete

# -------------------------------------------------- #
# 3. Prepare set of significant and background genes #
# -------------------------------------------------- #

# Extract unique contrasts
contrasts <- unique(deseq$contrast)

sign_genes_list <- list()

for (contrast in contrasts) {
  sign_genes <- deseq %>% 
    filter(padj < padj_cutoff & abs(log2FoldChange) >= log2_FC_cutoff) %>% 
    pull(Gene) %>% 
    unique()
  
  sign_genes_list[[contrast]] <- sign_genes
}

# Define background genes (all detected genes in the analysis)
background_genes <- deseq %>% 
  pull(Gene) %>% 
  unique()

# ----------------------------------- #
# 4. Set Organism Annotation Database #
# ----------------------------------- #
OrgDb_used <- org.Hs.eg.db ### delete
organism_used <- "hsa" ### delete

if (opts$annotation == "m") {
  OrgDb_used <- org.Mm.eg.db
  organism_used <- "mm"
} else if (opts$annotation == "h") {
  OrgDb_used <- org.Hs.eg.db
  organism_used <- "hsa"
} else {
  stop("Invalid annotation option provided. Use 'm' for mouse or 'h' for human.", call.=FALSE)
}


# Convert gene names to ENTREZID format using the dynamically selected OrgDb
background_ids <- bitr(background_genes,
                       fromType = "SYMBOL",
                       toType = "ENTREZID",
                       OrgDb = OrgDb_used)


# ----------------------------------------------- #
# 5. Perform KEGG Over-representation Analysis (ORA) #
# ----------------------------------------------- #
KEGG_results <- list()
for(i in names(sign_genes_list)){
  gene_symbols <- sign_genes_list[[i]] 
  
  entrez_ids <- bitr(gene_symbols, 
                     fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = OrgDb_used)
  
  # Define the organism parameter for KEGG enrichment based on annotation
  organism_used <- ifelse(opts$annotation == "m", "mmu", "hsa")
  
  kegg <- enrichKEGG(gene = entrez_ids$ENTREZID, 
                     organism = organism_used, 
                     pvalueCutoff = 0.05, 
                     pAdjustMethod = "BH", 
                     universe = background_ids$ENTREZID) 
  
  # prepare the dataframe from the results for making the figure  
  if (opts$annotation == "m") {
    kegg_df <- as.data.frame(kegg) %>%
      mutate(Description = gsub(" - Mus musculus \\(house mouse\\)", "", Description)) %>%
      filter(p.adjust < 0.05) %>% 
      arrange(p.adjust) %>% 
      mutate(Description = factor(Description, levels = rev(unique(Description))),
             padj = ifelse(`p.adjust` < 0.05, "padj < 0.05", "padj > 0.05"),
             `-log10(padj)` = -log10(p.adjust))
  } else {
    kegg_df <- as.data.frame(kegg) %>%
      filter(p.adjust < 0.05) %>% 
      arrange(p.adjust) %>% 
      mutate(Description = factor(Description, levels = rev(unique(Description))),
             padj = ifelse(`p.adjust` < 0.05, "padj < 0.05", "padj > 0.05"),
             `-log10(padj)` = -log10(p.adjust))
  }
  
  # Save the kegg_df results as an Excel file
  df_file <- file.path(opts$output, paste0("KEGG_", contrast, "_results.xlsx"))
  writexl::write_xlsx(kegg_df, df_file)
  
  p_kegg <- ggplot(head(kegg_df, 15), aes(x = `-log10(padj)`, y = Description)) +
    geom_bar(stat = "identity", fill = "grey90", color = "black") +
    geom_text(aes(label = Count), vjust = 0.5, hjust = -0.5) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed") +
    scale_x_continuous(limits = c(0, max(kegg_df$`-log10(padj)` + .5))) +
    labs(x = "-log10(padj)",
         y = "Enriched Term",
         size = "#Genes",
         fill = "") +
    ggtitle(paste0("KEGG: ", i)) + 
    my_theme + 
    theme(legend.title = element_text(),
          legend.position = c(0.7, 0.15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  
  # dynamic height based on number of enriched pathways
  n_pathways <- length(unique(kegg_df$Description))
  base_height <- 1.3
  height_pkegg <- if (n_pathways == 1) base_height else base_height + n_pathways * 0.2
  
  output_file_pkegg <- file.path(output_dir, paste0("KEGG_", contrast, ".pdf"))
  save_ggplot(p_kegg, output_file_pkegg, width = 5, height = height_pkegg)
  message("Saved KEGG plot to: ", output_file_pkegg)
  
}

message("done with KEGG analysis")
 





# --------------------------------- #
# 6. Perform GO Enrichment Analysis #
# --------------------------------- #
# Define the GO enrichment function (produces both an Excel file and a plot)
Go_plot <- function(genes, gene_map, background, ont, output_dir, contrast) {
  go <- enrichGO(gene = genes,
                 universe = background,
                 OrgDb = OrgDb_used,
                 ont = ont,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH")
  
  # Check if any enrichment results were produced
  if (is.null(go) || nrow(go@result) == 0) {
    message("No GO enrichment results for ", contrast, " in ", ont, " ontology.")
    return(NULL)
  }
  
  # Create an expanded results dataframe (one row per gene)
  go_df <- go@result %>% 
    as.data.frame() %>% 
    filter(p.adjust < 0.05) %>% 
    dplyr::rename(ENTREZID = geneID) %>% 
    separate_rows(ENTREZID, sep = "/") %>% 
    left_join(gene_map, by = "ENTREZID") %>% 
    dplyr::rename(Gene = SYMBOL) %>% 
    arrange(p.adjust)
  
  # Save the GO results dataframe as an Excel file
  excel_filename <- file.path(output_dir, paste0("GO_", ont, "_", contrast, "_results.xlsx"))
  write_xlsx(go_df, excel_filename)
  message("Saved GO enrichment results (Excel) to: ", excel_filename)
  
  # Select the top 10 pathways (by adjusted P-value) for the plot
  top_pathways <- go@result %>% 
    filter(p.adjust < 0.05) %>% 
    group_by(Description) %>% 
    summarise(min_padj = min(p.adjust)) %>% 
    arrange(min_padj) %>% 
    pull(Description) %>% 
    unique() %>% 
    head(10)
  
  top10 <- go@result %>% filter(Description %in% top_pathways)
  
  # Reorder the factor levels for plot display
  top10$Description <- factor(top10$Description, levels = top10$Description)
  
  # Create the horizontal bar plot
  p_go <- ggplot(top10, aes(x = p.adjust, y = reorder(Description, -p.adjust))) +
    geom_bar(stat = "identity", color = "black", fill = "grey80") +
    geom_text(aes(label = Count), vjust = 0.5, hjust = -0.5) +
    labs(x = "Adjusted P-Value", y = "") +
    ggtitle(if (ont == "BP") {
      "GO Biological Process"
    } else if (ont == "MF") {
      "GO Molecular Function"
    } else {
      "GO Cellular Component"
    }) +
    geom_vline(xintercept = 0.05, linetype = "dashed") +
    scale_x_continuous(trans = trans_reverser("log10")) +
    my_theme
  
  # Adjust the plot height dynamically based on the number of enriched pathways
  n_pathways <- length(unique(top10$Description))
  base_height <- 1.3
  height_p_go <- if (n_pathways == 1) base_height else base_height + n_pathways * 0.2
  
  # Define the output file name for the PDF plot
  pdf_filename <- file.path(output_dir, paste0("GO_", ont, "_", contrast, ".pdf"))
  save_ggplot(p_go, pdf_filename, width = 6.5, height = height_p_go)
  message("Saved GO plot to: ", pdf_filename)
  
  return(go_df)
}


# Loop over each contrast (gene list) and each ontology (BP, MF, CC)
message("Performing GO enrichment analysis")
for (contrast in names(sign_genes_list)) {
  gene_symbols <- sign_genes_list[[contrast]]
  
  # Convert the significant genes for this contrast into ENTREZ IDs and create a mapping table
  gene_map <- bitr(gene_symbols, 
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = OrgDb_used)
  
  if (nrow(gene_map) == 0) {
    message("No valid ENTREZ IDs found for contrast: ", contrast)
    next
  }
  
  for (ont in c("BP", "MF", "CC")) {
    message("Running GO enrichment for contrast: ", contrast, " | Ontology: ", ont)
    Go_plot(genes = gene_map$ENTREZID,
            gene_map = gene_map,
            background = background_ids$ENTREZID,
            ont = ont,
            output_dir = output_dir,
            contrast = contrast)
  }
}

message("Done with GO enrichment analysis")


 
# ############################################################
# ######## Perform Wikipathways Enrichment Analysis ##########
# ############################################################
# wiki <- enrichWP(gene=entrez_ids$ENTREZID,
#                  organism = "Mus musculus") # set the organism
# 
# 
# # Alternatively (if the wikipathway data cannot be retrieved automatically)
# download.file("https://data.wikipathways.org/current/gmt/wikipathways-20240810-gmt-Mus_musculus.gmt", destfile = "Mus_musculus.gmt") # check website and update to the latest version
# 
# # Load the downloaded GMT file
# gmt_data <- read.gmt("Mus_musculus.gmt")
# 
# # Enrich using the GMT data
# wiki <- enricher(gene = entrez_ids$ENTREZID, TERM2GENE = gmt_data)
# 
# 
# # transform the wiki object to a dataframe, sort by adjusted pvalue and create a new column "padj"
# wiki_df <- as.data.frame(wiki) %>%
#   arrange(p.adjust) %>% 
#   mutate(padj = ifelse(`p.adjust` < 0.05, "padj < 0.05", "padj > 0.05"),
#          `-log10(padj)` = -log10(p.adjust))
# 
# 
# # adjust the Description column (when organism = "Mus musculus")
# wiki_df$Description = sub("%.*", "", wiki_df$Description)
# 
# # set the Description column as a factor in reverse Order (so the pathways are ranked from top to bottom in the ggplot figure)
# wiki_df$Description = factor(wiki_df$Description, levels = rev(unique(wiki_df$Description)))
# 
# # function to generate the plot
# wiki_plot <- function(df){
#   ggplot(head(df, 25), aes(x = `-log10(padj)`, y = Description)) + # adapt the head() function if you want to show other than 25 pathways
#     geom_bar(stat = "identity", fill = "grey90", color = "black") +
#     geom_text(aes(label = Count), vjust = 0.5, hjust=-0.5) +
#     geom_vline(xintercept = -log10(0.05), linetype = "dashed") +
#     scale_x_continuous(limits = c(0, max(df$`-log10(padj)`+.5))) +
#     labs(x = "-log10(padj)",
#          y = "Enriched Term",
#          fill = "") +
#     ggtitle("Wikipathways") +
#     my_theme +
#     theme(legend.title = element_text(),
#           legend.position = c(0.7,0.15)) # adapt accordingly (or remove with legend.position = "None")
# }
# 
# 
# p <- wiki_plot(wiki_df)
# 
# 
# #save the plot as pdf file
# plot_name <- "[date]_[initials]_[project]_[whatever you want].pdf"
# ggsave(plot_name, plot = p, width = 7, height = 3) # adapt width and height accordingly
# 
# 
# 
# 
# 
# 
# 

# 
# 
# ######################################################
# ######## Gene-set enrichment analysis (GSEA) #########
# ######################################################
# 
# # create the df sign_genes containing the genes of interest in decreasing order (here: DESeq2 result dataframe sorted by log2FoldChange)
# genes <- df$log2FoldChange # df = dataframe containing the results
# names(genes) <- df$geneID
# genes <- sort(genes, decreasing = TRUE)
# 
# 
# Go_gsea_plot <- function(input_genes, ont, plot = TRUE) {
#   # Remove duplicate gene names
#   input_genes <- input_genes[!duplicated(names(input_genes))]
#   
#   # Sort the genes in decreasing order
#   input_genes <- sort(input_genes, decreasing = TRUE)
#   
#   # Run gseGO with verbose output
#   go <- gseGO(
#     geneList = input_genes,
#     ont = ont,
#     keyType = "ENSEMBL",
#     OrgDb = org.Mm.eg.db,
#     pvalueCutoff = 0.05,
#     pAdjustMethod = "BH",
#     minGSSize = 5,
#     verbose = TRUE
#   )
#   
#   # Check if any gene sets are enriched
#   if (is.null(go) || nrow(go@result) == 0) {
#     message("No significant gene sets were found.")
#     return(list(result = NULL, leading_edge_genes = NULL))
#   }
#   
#   result <- go@result
#   # leading_edge_genes <- leading_edge(go, geneSetID = 1)
#   
#   # Plotting the results using gseaplot2 from enrichplot package if plot = TRUE
#   if (plot) {
#     gseaplot2(go, geneSetID = 1)
#   }
#   
#   # Return the result list with the leading edge genes
#   return(list(result = result)) #, leading_edge_genes = leading_edge_genes))
# }
# 
# 
# ontologies <- c("BP", "MF", "CC")
# p_results <- list() # create an empty list to store the GSEA results for each ontology
# 
# # run a for-loop to do the enrichment for each ontology  
# for(ont in ontologies){
#   p <- Go_gsea_plot(genes, ont, plot = FALSE)
#   
#   p_results[[ont]] <- p
# }
# 
# 
# # Loop through the list and save each plot using ggsave
# for (i in seq_along(p_results)) {
#   plot_name <- paste0("GeneSet_enrichment_", names(p_results)[i], ".pdf")
#   ggsave(plot_name, plot = p_results[[i]], width = 8.5, height = 6) # adapt accordingly
# }