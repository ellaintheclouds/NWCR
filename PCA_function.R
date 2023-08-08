#0 Set up-----------------------------------------------------------------------

setwd("C:/Users/richarej/tcga//app versions/tcga_app_2_07.08.23")

## Loading libraries ##
if (!require("ggcorrplot", quietly = TRUE)) { install.packages("ggcorrplot") } 
if (!require("factoextra", quietly = TRUE)) { install.packages("factoextra") } 
if (!require("ggfortify", quietly = TRUE)) { install.packages("ggfortify") } 
if (!require("umap", quietly = TRUE)) { install.packages("umap") } 
if (!require("plotly", quietly = TRUE)) { install.packages("plotly") } 
if (!require("viridis", quietly = TRUE)) { install.packages("viridis") }

# General
library(stringr) # Makes working with strings simpler.
library(viridis) # Colour scheme.

# Biological data
library(TCGAbiolinks) # For integrative analysis of GDC data.
library(ggcorrplot) # Allows visulaisation of a correlation matrix using ggplot2.
library(SummarizedExperiment) # For storing processed data from high-throughput sequencing assays.
library(edgeR) # Differential expression analysis of RNA-seq expression profiles.
# Allows reordering  of matrix correlation and displays significance level on the plot. Allows computation of a matrix of p-values.
library(plotly) # Creates interactive, publication-suitable graphs.

# PCA
library(factoextra) # For principal component analysis
library(ggfortify) # Allows plotting of PCA and survival analysis.
library(umap) # Algorithm for dimensional reduction


#1 Start of function------------------------------------------------------------
pca_function <- function(cancer_type_list, gene_string){
  gene_list <- strsplit(gene_string, split = "\n")[[1]] # Split the user input into a list

  # storing all plots
  out_plots <- list()
  
#2 Formatting and looping through cancer type data------------------------------
   for(current_cancer_type in cancer_type_list){
    
    current_file_name <- paste0(current_cancer_type, ".rds")
    print(current_file_name) # Prints the names of the files that are being accessed in the loop
    
    out_plots[[current_cancer_type]] <- list()
    
    #data_RNAseq <- readRDS(current_file_name) # Reads in the file to be used in this iteration of the loop
    data_RNAseq <- readRDS(paste0("data/", current_file_name))  

  #3 Get the count data and the clinical annotation-------------------------------
    count_data <- assay(data_RNAseq, "unstranded") # gene count matrix
    tpm_matrix <- assay(data_RNAseq, "tpm_unstrand") # transcript per million matrix (for plotting colour later)
    
    sample_annotation <- colData(data_RNAseq) # clinical annotation data
    gene_annotation <- rowData(data_RNAseq) # gene annotation data
    
    # Normalise data
    x <- DGEList(counts = count_data, genes = gene_annotation)
    keep.exprs <- filterByExpr(x)
    x <- x[keep.exprs,, keep.lib.sizes=FALSE]
    x <- calcNormFactors(x, method = "TMM")
    
    lcpm_filtered <- cpm(x, log=TRUE)
    
    
    #4 PCA------------------------------------------------------------------------
    #gene_annotation <- read.csv("/Count and annotation data/gene_annotation.csv") # CHANGE to add gene annotations
    #count_data <- read.csv("/Count and annotation data/count_data.csv", row.names=1) # CHANGE to create count matrix
    #count_data <- data.matrix(count_data) # converted to a matrix as these are good for calculations
    count_pca <- prcomp(t(lcpm_filtered), scale=TRUE) # Carry out PCA on transformed matrix of count data.
    fviz_eig(count_pca, addlabels = TRUE) # Visualise the eiganvalues.
    
    
    #5 Reorganise PCA output for plotting-----------------------------------------
    
    # Reorganise to calculate the strength of contribution.
    gene_contribution <- as.data.frame(count_pca$rotation)# (number of PC equal to number of variables in the input)
    
    gene_contribution_abs <- abs(gene_contribution) # abs() (absolute values) makes every value positive
    
    gene_contribution_rank <- apply(gene_contribution_abs, 2, rank) # rank the contribution
    gene_contribution_percentile <- as.data.frame(gene_contribution_rank/nrow(gene_contribution_rank)) # create percentiles, with the lower values indicating higher contribution
    
    # Putting column and row names back in after having been in matrix form.
    pc_names <- colnames(count_pca$rotation)
    gene_names <- rownames(count_pca$rotation)
    colnames(gene_contribution) <- pc_names
    rownames(gene_contribution) <- gene_names
    colnames(gene_contribution_percentile) <- pc_names
    rownames(gene_contribution_percentile) <- gene_names
    
    # The x-/y-coordinates to plot for each PC
    pca_positions <- as.data.frame(count_pca$x)
    pca_summary <- summary(count_pca)$importance[2,] * 100 # Showing the contribution of each position, to represent its importance.
    

    #6 Formatting and looping through gene data-------------------------------------
    for(current_gene in gene_list){
      
      out_plots[[current_cancer_type]][[current_gene]] <- list()
      
      tpm_dataframe <- as.data.frame(t(tpm_matrix))
      tpm_dataframe$sample_id_tpm <- rownames(tpm_dataframe)
      tpm_dataframe <- tpm_dataframe[, c("sample_id_tpm", current_gene)]
      colnames(tpm_dataframe) <- c("sample_id_tpm", "tpm")
      
      sample_annotation_tpm <- merge(sample_annotation, tpm_dataframe, by = 0)
      row.names(sample_annotation_tpm) <- sample_annotation_tpm$sample_id_tpm
      sample_annotation_tpm$Row.names <- NULL
    
      # Gene relationships
      # pca_plot_df <- as.data.frame(merge(gene_annotation, as.data.frame(count_pca$x), by.x= 0, by.y= 0))
      pca_plot_df <- as.data.frame(merge(sample_annotation_tpm, as.data.frame(count_pca$x), by.x= 0, by.y= 0))
        
        
      #7 Start plot function for every possible combination of dimensions-------------
        for(pcx in 1:1){
          for(pcy in 1:2){
            
            
      #8 Plotting PCA-----------------------------------------------------------------
            # the if statement is to remove duplicate plots (i.e. it will just plot PC3 vs PC4, and not PC4 vs PC3)
            if(pcx < pcy){
              plotting_dimentions <- c(pcx, pcy)
              
              ## Find top contributing genes for pcx
              top_contrib_genes_x <- data.frame(contribution = gene_contribution[,pcx], 
                                                gene_contribution_percentile = gene_contribution_percentile[,pcx], 
                                                gene = row.names(gene_contribution))
              # order by contribution percentile
              top_contrib_genes_x <- top_contrib_genes_x[order(top_contrib_genes_x$gene_contribution_percentile),]
              print(paste0("Top contributing genes for ", pcx, ":"))
              print(head(top_contrib_genes_x))
              
              ## Find top contributing genes for pcy
              top_contrib_genes_y <- data.frame(contribution = gene_contribution[,pcy], 
                                                gene_contribution_percentile = gene_contribution_percentile[,pcy], 
                                                gene = row.names(gene_contribution))
              # order by contribution percentile
              top_contrib_genes_y <- top_contrib_genes_y[order(top_contrib_genes_y$gene_contribution_percentile),]
              print(paste0("Top contributing gene for ", pcy, ":"))
              print(head(top_contrib_genes_y))
              
              ## plot PCA
              current_pca_plot_df <- pca_plot_df[,c("sample_submitter_id", "tpm", paste0("PC", pcx), paste0("PC", pcy))]
              colnames(current_pca_plot_df) <- c("sample_submitter_id", "tpm", "pcx", "pcy")
              
              # Labels for the axis
              tpm_label <- "\n Expression \n (TPM)"
              
              # Plot
              pca_plot <- ggplot(data = current_pca_plot_df, aes(x = pcx, y = pcy)) + 
                geom_point(aes(colour = tpm)) +
                #scale_x_continuous(limits = c(-160, 100), breaks = c(-200, -150, -100, -50, 0, 50, 100)) +
                #scale_y_continuous(limits = c(-40, 40), breaks = c(-40, -30, -20, -10, 0, 10, 20, 30, 40)) + 
                labs(colour = paste0(current_gene, tpm_label)) + 
                theme_bw() + 
                xlab(paste0("PC", pcx, " (", pca_summary[pcx], ")")) + 
                ylab(paste0("PC", pcy, " (", pca_summary[pcy], ")")) + 
                scale_colour_viridis()
              
              out_plots[[current_cancer_type]][[current_gene]][[paste0(pcx, "_", pcy)]] <- pca_plot
              
              #dir.create(file.path(paste0("out/pca/", current_cancer_type)), showWarnings = FALSE)
              #dir.create(file.path(paste0("out/pca/", current_cancer_type, "/", current_gene)), showWarnings = FALSE)
              
              #ggsave(paste0("out/pca/", current_cancer_type, "/", current_gene, "/PC", pcx, "_PC", pcy, ".png"), pca_plot, height = 12, width = 16, units = "cm")              
            } # Plotting PCs when pcx < pcy
          } # pcy
        } # pcx
    } # Gene
   } # Cancer type
  
  return(out_plots)
  
} # Function

#pca_function(c("TCGA-LAML", "TCGA-ACC"), "ENSG00000000003.15")