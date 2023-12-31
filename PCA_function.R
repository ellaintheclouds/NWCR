#0 Set up----------------------------------------------------------------------
# Turn this to true if loading outside the shinyapp
if(FALSE){
  source("script/load_libraries.R")
  source("script/load_variables.R")
}


#1 Start of function------------------------------------------------------------
pca_function <- function(cancer_type_list, gene_list){
  
  formatting_output <- input_format(gene_list)
  
  # storing all plots
  output_data <- list("pca plots", "contribution plots", "contribution percentile dataframes", "pca dataframes")
  output_data[["contribution plots"]] <- list()
  output_data[["pca plots"]] <- list()
  
  
  #2 Formatting and looping through cancer type data------------------------------
  for(current_cancer_type in cancer_type_list){
    
    # Re-assigning names to cancer types
    current_cancer_type_named <- current_cancer_type 
    names(current_cancer_type) <- list_of_cancer_types_rev[current_cancer_type]
    
    current_file_name <- paste0(current_cancer_type, ".rds")
    print(current_file_name) # Prints the names of the files that are being accessed in the loop
    
    output_data[["pca plots"]][[current_cancer_type]] <- list()
    
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
    
    # Making scree plots----------
    # Foratting data: shorten, transform, merge
    pca_summary <- summary(count_pca)$importance[2,] * 100 # Showing the contribution of each position, to represent its importance.----------------------------------------------------------------------------------------------------------------------------------------------------
    contribution_dataframe <- data.frame(PC = paste0("PC", 1:10), contribution = pca_summary[1:10])
    percentile_dataframe <- as.data.frame(t(as.data.frame(gene_contribution_percentile)))
    percentile_dataframe$PC <- paste0("PC", 1:nrow(percentile_dataframe))
    
    contr_and_per_df <- merge(contribution_dataframe, percentile_dataframe, by = "PC")
    contr_and_per_df$PC <- factor(contr_and_per_df$PC, levels = paste0("PC", 1:10))
    contr_and_per_df <- contr_and_per_df[order(contr_and_per_df$PC),]
    row.names(contr_and_per_df) <- contr_and_per_df$PC
    
    output_data[["contribution percentile dataframes"]][[current_cancer_type]] <- contr_and_per_df
    
    
    #6 Formatting and looping through gene data-------------------------------------
    for(current_gene in gene_list){
      
      output_data[["pca plots"]][[current_cancer_type]][[current_gene]] <- list()
      
      # Current gene name
      current_gene_name <- formatting_output[["formatted_gene_list"]][formatting_output[["formatted_gene_list"]]$gene_id == current_gene, "gene_name"]
      
      # Tpm 
      tpm_dataframe <- as.data.frame(t(tpm_matrix))
      tpm_dataframe$sample_id_tpm <- rownames(tpm_dataframe)
      tpm_dataframe <- tpm_dataframe[, c("sample_id_tpm", current_gene)]
      colnames(tpm_dataframe) <- c("sample_id_tpm", "tpm")
      
      sample_annotation_tpm <- merge(sample_annotation, tpm_dataframe, by = 0)
      row.names(sample_annotation_tpm) <- sample_annotation_tpm$sample_id_tpm
      sample_annotation_tpm$Row.names <- NULL
      
      # Gene relationships
      pca_plot_df <- as.data.frame(merge(sample_annotation_tpm, as.data.frame(count_pca$x), by.x= 0, by.y= 0))
      
      
      #7 Start plot function for every possible combination of dimensions-------------
      for(pcx in 1:10){
        for(pcy in 1:10){
          
          
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
            tpm_label <- "Transcripts\nper\nMillion"
            
            pca_title <- paste0(current_gene_name, " Expression in ", names(current_cancer_type))
            
            # Plot
            pca_plot <- ggplot(data = current_pca_plot_df, aes(x = pcx, y = pcy)) + 
              geom_point(aes(colour = tpm)) +
              labs(title = pca_title, colour = tpm_label) + 
              theme_minimal() + 
              xlab(paste0("PC", pcx, " (", pca_summary[pcx], ")")) + 
              ylab(paste0("PC", pcy, " (", pca_summary[pcy], ")")) + 
              scale_color_viridis(option = "A")
            
            output_data[["pca plots"]][[current_cancer_type]][[current_gene]][[paste0(pcx, "_", pcy)]] <- pca_plot
          } # Plotting PCs when pcx < pcy
        } # pcy
      } # pcx
    } # Gene
  } # Cancer type
  
  return(output_data)
  
} # Function
