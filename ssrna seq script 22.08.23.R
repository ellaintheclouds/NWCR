### MERKEL CELL ###
# Process the merkel cells from cellranger output for the fingerprint study

# Single nuclei data (snRNA seq) 
# 10x genomics dataset (droplet-based sequencing)
# Amplifies polyA (only studies spliced RNA in practice, but in reality, it has some non-specific binding)

# Set up------------------------------------------------------------------------
setwd("E:/NWCR/R working directory/Week 8/")

# Installation
if(FALSE){
  BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                         'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                         'SummarizedExperiment', 'batchelor', 'HDF5Array',
                         'terra', 'ggrastr'))
  
  install.packages("devtools")
  install.packages('terra', repos='https://rspatial.r-universe.dev')
  install.packages('R.utils')
  install.packges("dplyr")
  
  devtools::install_github('cole-trapnell-lab/monocle3')
  BiocManager::install("celda")
  BiocManager::install("scater")
  BiocManager::install("celldex")
  BiocManager::install("scDblFinder")
} # If false

# Import libraries
library(monocle3)
library(ggplot2)
library(scDblFinder)
library(scater)
library(celldex)
library(magrittr)
library(dplyr)


# Get data----------------------------------------------------------------------
cellranger_output <- c("Data/ventral 1", "Data/ventral 2")

cds <- lapply(cellranger_output, function(x)load_mm_data(mat_path = paste0(x, "/matrix.mtx"), 
                                                         feature_anno_path = paste0(x, "/features.tsv"), 
                                                         cell_anno_path = paste0(x, "/barcodes.tsv")))

colnames(rowData(cds[[1]])) <- c("gene_short_name", "gene_type")
colnames(rowData(cds[[2]])) <- c("gene_short_name", "gene_type")

# Provide names for data
names(cds) <- gsub("Data/", "", cellranger_output)


# QC filter each sample individually----------
for(idx in 1:length(cds)){
  
  # Finding mitochondrial and ribosomal genes
  is_mito <- grepl("^MT",rowData(cds[[idx]])$gene_short_name) #^ means starts with and [] means either or
  is_ribo <- grepl("^RP[SL]",rowData(cds[[idx]])$gene_short_name) 
  
  # Compute per-cell quality control metrics for a count matrix/single cell experiment
  umi_cell <- perCellQCMetrics(cds[[idx]],subsets=list(Mito=is_mito, Ribo = is_ribo))
  
  # Adding names to cds
  for (name in colnames(umi_cell)){
    cds[[idx]][[name]] = umi_cell[[name]]
  } # For name
  cds[[idx]]$sample_id <- names(cds)[idx]
} # For idx

head(colData(cds[[1]]))

# Creating a data frame
metadata <- as.data.frame(do.call(rbind, lapply(cds, colData)))

# Setting bounds----------
abline_val_lower1 <- c(2000, 2000, 1200, 0, 0, 0, 2000)
abline_val_lower2 <- c(2000, 2000, 1200, 0, 0, 0, 2000)
abline_val_upper1 <- c(50000, 50000, 8000, 0, 0, 0, 30000)
abline_val_upper2 <- c(50000, 50000, 8000, 0, 0, 0, 30000)

plot_val <- c("Number of UMI" = "n.umi", "Sum" = "sum", "Genes Detected" = "detected", 
              "Mitochondrial Subset Sum" = "subsets_Mito_sum", 
              "Mitochondrial Subset Genes Detected" = "subsets_Mito_detected", 
              "Mitochondrial Subset Percent" = "subsets_Mito_percent", "Total" = "total")

# Look at overall statistics----------
for(idx in 1:length(plot_val)){
  val <- plot_val[idx]
  
  # Creating a plot
  p <- eval(parse(text = paste0("ggplot(metadata, aes( x = ", val, "))"))) + 
    geom_density(alpha = 0.2, aes(color=sample_id, fill= sample_id)) + 
    scale_x_log10() + 
    geom_vline(xintercept = abline_val_lower1[idx], colour="pink") +
    geom_vline(xintercept = abline_val_upper1[idx], colour="pink") +
    geom_vline(xintercept = abline_val_lower2[idx], colour="deepskyblue2") +
    geom_vline(xintercept = abline_val_upper2[idx], colour="deepskyblue2") +
    xlab("Total") + ylab("Density") + 
    theme_classic() +
    ggtitle(names(val))
  
  print(p)
  ggsave(filename = paste0(names(val), ".png"), 
         plot = p, 
         path = "Output/Metadata", 
         width = 8, height = 6)
} # For idx

# Remove cells that failed quality control parameters----------
# (manually decided on discard threshold by viewing the graph)
thresholds_lower <- list(n.umi=c('ventral 1'=1000, 'ventral 2'=2000), 
                         detected = c('ventral 1'=1000, 'ventral 2'=1200))

thresholds_upper <- list(n.umi=c('ventral 1'=3000, 'ventral 2'=5000), 
                         detected = c('ventral 1'=10000, 'ventral 2'=10000))

# Remove cells with low counts or low umi based on i
cds_filtered <- list()
for(idx in 1:length(cds)){
  
  # Get the associated values for the current sample
  current_sample <- names(cds)[idx]
  current_threshold_lower_n_umi <- thresholds_lower[["n.umi"]][[current_sample]]
  current_threshold_upper_n_umi <- thresholds_upper[["n.umi"]][[current_sample]]
  current_threshold_lower_detected <- thresholds_lower[["detected"]][[current_sample]]
  current_threshold_upper_detected <- thresholds_upper[["detected"]][[current_sample]]
  
  current_cds <- cds[[current_sample]]
  
  # Subset data removing cells with low n.umi or low genes detected
  cds_filtered[[current_sample]] <- 
    current_cds[ ,	  # If necessary, rowData filter goes here before the comma
                 colData(current_cds) %>%
                   subset(
                     (n.umi >  current_threshold_lower_n_umi) & 
                       (n.umi <  current_threshold_upper_n_umi) & 
                       (detected >  current_threshold_lower_detected) &
                       (detected <  current_threshold_upper_detected)
                   ) %>%
                   row.names
    ] # current_cds
}# For

names(cds_filtered) <- c("ventral 1", "ventral 2")

# combine all samples----------
cds_all <- combine_cds(cds_filtered, sample_col_name ="sample_id")


# Preprocessing-----------------------------------------------------------------
set.seed(123) # Used to create reproducible results by creating variables that take on random valuables.

# preprocess the data, dimention redution and clustering
cds_all <- preprocess_cds(cds_all, num_dim = 100)
cds_all <- reduce_dimension(cds_all)
cds_all <- cluster_cells(cds_all)

# find doublets with scDblFinder, use the clusters identified for doublet finding
set.seed(123)
cds_all$clusters <- clusters(cds_all)

# find clusters that are doublets
# see https://rdrr.io/github/plger/scDblFinder/man/findDoubletClusters.html for explaination on the outputs
tab <- findDoubletClusters(cds_all, cds_all$clusters)# Doublet cluster

# classify doublet "cells"
sce <- scDblFinder(cds_all, samples = "sample_id", dbr = 0.1, clusters = "clusters") # Doublet label
cds_all$scDblFinder <- sce$scDblFinder.class

# How to specify the number of dimensions/regress unwanted variation
# https://github.com/cole-trapnell-lab/monocle-release/issues/178

# Plot reduced------------------------------------------------------------------

# Reduce
cds_all <- reduce_dimension(cds_all)

# Plot reduced
cluster_umap <- plot_cells(cds_all)
ggsave(filename = "Reduced.png", plot = cluster_umap, path = "Output/UMAP/Reduced", 
                          width = 8, height = 6)

red_cluster_umap <- plot_cells(cds_all, color_cells_by = "clusters")
ggsave(filename = "Reduced Cluster.png", plot = red_cluster_umap, path = "Output/UMAP/Reduced", 
       width = 8, height = 6)

# Plot by genes
red_cluster_KRT20_KRT10_KRT1_umap <- plot_cells(cds_all, color_cells_by = "clusters", genes = c("KRT20", "KRT10", "KRT1"))
ggsave(filename = "Reduced Cluster KRT20 KRT10 KRT1.png", plot = red_cluster_KRT20_KRT10_KRT1_umap, path = "Output/UMAP/Reduced", 
       width = 8, height = 6)

red_cluster_COL1A1_TNC_umap <- plot_cells(cds_all, color_cells_by = "clusters", genes = c("COL1A1", "TNC"))
ggsave(filename = "Reduced Cluster COL1A1 TNC.png", plot = red_cluster_COL1A1_TNC_umap, path = "Output/UMAP/Reduced", 
       width = 8, height = 6)

red_cluster_HBB_umap <- plot_cells(cds_all, color_cells_by = "clusters", genes = "HBB")
ggsave(filename = "Reduced Cluster HBB.png", plot = red_cluster_HBB_umap, path = "Output/UMAP/Reduced", 
       width = 8, height = 6)

# Plot by doublet
red_doublet_umap <- plot_cells(cds_all, color_cells_by = "scDblFinder")
ggsave(filename = "Reduced Doublet.png", plot = red_doublet_umap, path = "Output/UMAP/Reduced", 
       width = 8, height = 6)

# Reduce again
cds_all2 <- reduce_dimension(cds_all, reduction_method = "tSNE")

# Plot double reduced 
red2_umap <- plot_cells(cds_all2, reduction_method = "tSNE", colour_cells_by = "clusters")
ggsave(filename = "Reduced2.png", plot = red2_umap, path = "Output/UMAP/Reduced", 
       width = 8, height = 6)


# Plot with order---------------------------------------------------------------
# Order
cds_all <- learn_graph(cds_all)
cds_all <- order_cells(cds_all)

ord_cluster_umap <- plot_cells(cds_all, color_cells_by = "clusters")
ggsave(filename = "Ordered Cluster.png", plot = ord_cluster_umap, path = "Output/UMAP/Ordered", 
       width = 8, height = 6)

# Plot by genes
ord_cluster_KRT20_KRT10_KRT1_umap <- plot_cells(cds_all, color_cells_by = "clusters", genes = c("KRT20", "KRT10", "KRT1"))
ggsave(filename = "Ordered Cluster KRT20 KRT10 KRT1.png", plot = ord_cluster_KRT20_KRT10_KRT1_umap, path = "Output/UMAP/Ordered", 
       width = 8, height = 6)

ord_cluster_COL1A1_TNC_umap <- plot_cells(cds_all, color_cells_by = "clusters", genes = c("COL1A1", "TNC"))
ggsave(filename = "Ordered Cluster COL1A1 TNC.png", plot = ord_cluster_COL1A1_TNC_umap, path = "Output/UMAP/Ordered", 
       width = 8, height = 6)

ord_cluster_HBB_umap <- plot_cells(cds_all, color_cells_by = "clusters", genes = "HBB")
ggsave(filename = "Ordered Cluster HBB.png", plot = ord_cluster_HBB_umap, path = "Output/UMAP/Ordered", 
       width = 8, height = 6)

# Plot by doublet
ord_doublet_umap <- plot_cells(cds_all, color_cells_by = "scDblFinder")
ggsave(filename = "Reduced Doublet.png", plot = ord_doublet_umap, path = "Output/UMAP/Ordered", 
       width = 8, height = 6)

# Genes-------------------------------------------------------------------------
#marker_test_res <- top_markers(cds_all, group_cells_by = "cluster", 
#                               reference_cells = 1000, cores = 8)
#
#top_specific_markers <- marker_test_res %>%
#  filter(fraction_expressing >= 0.10) %>%
#  group_by(cell_group) %>%
#  top_n(1, pseudo_R2)
#
#top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))
#
#plot_genes_by_group(cds_all,
#                    top_specific_marker_ids,
#                    group_cells_by="partition",
#                    ordering_type="maximal_on_diag",
#                    max.size=3)
#
#top_specific_markers <- marker_test_res %>%
#  filter(fraction_expressing >= 0.10) %>%
#  group_by(cell_group) %>%
#  top_n(3, pseudo_R2)
#
#top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))
#
#plot_genes_by_group(cds_all,
#                    top_specific_marker_ids,
#                    group_cells_by="partition",
#                    ordering_type="cluster_row_col",
#                    max.size=3)
