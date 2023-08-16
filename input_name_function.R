# 0 Set up----------------------------------------------------------------------
setwd("D:/app versions/tcga_app_usb")

gene_annotation <- read.csv("data/gene_annotation.csv")


# 1 Function--------------------------------------------------------------------
input_format <- function(gene_list){
  formatting_output <- list("rejected_gene_list" = vector(), "formatted_gene_list" =list())
  idx <- 1
  for (g in gene_list){
    
    if (g %in% gene_annotation$gene_id | g %in% gene_annotation$gene_id_stable |g %in%  gene_annotation$gene_name){
      
      gene_annotation_relavent <- gene_annotation[gene_annotation$gene_id == g | gene_annotation$gene_id_stable == g | gene_annotation$gene_name ==g, c("gene_id", "gene_name")]
      gene_annotation_relavent$merged_name <- paste0(gene_annotation_relavent$gene_name, ": ", gene_annotation_relavent$gene_id)
      formatting_output[["formatted_gene_list"]][[idx]] <-   gene_annotation_relavent
      
      idx <- idx+1
    } # If
    
    else {
      formatting_output[["rejected_gene_list"]] <- c(formatting_output[["rejected_gene_list"]], g)
    } # Else
  } # For
  
  formatting_output[["formatted_gene_list"]] <- do.call(rbind, formatting_output[["formatted_gene_list"]])
  return(formatting_output)
  
} # function
