# 0 Set up----------------------------------------------------------------------
setwd("D:/app versions/tcga_app_usb")

gene_annotation <- read.csv("data/gene_annotation.csv")

# 1 Function start--------------------------------------------------------------
input_format <- function(gene_list){
  formatting_output <- list("rejected_gene_list" = vector(), "formatted_gene_list" =vector(), "gene_names_assigned" = vector())
  
  for (g in gene_list){
    
    if (g %in% gene_annotation$gene_id | g %in% gene_annotation$gene_id_stable |g %in%  gene_annotation$gene_name){
      
      gene_annotation_relavent <- gene_annotation$gene_id[gene_annotation$gene_id == g | gene_annotation$gene_id_stable == g | gene_annotation$gene_name ==g]
      formatting_output[["formatted_gene_list"]] <- c(formatting_output[["formatted_gene_list"]], gene_annotation_relavent)
      
      gene_name_assigned <- gene_annotation$gene_name[gene_annotation$gene_id == g | gene_annotation$gene_id_stable == g | gene_annotation$gene_name ==g]
      formatting_output[["gene_names_assigned"]] <- c(formatting_output[["gene_names_assigned"]], gene_name_assigned)
      
    } # if in gene_annotation
    
    else {
      formatting_output[["rejected_gene_list"]] <- g
    } # Else
  } # For
  
  return(formatting_output)
  
} # function

#formatting_output <- input_format(c("ENSG00000000419.13", 
#"ENSG00000000457.14", 
#"ENSG00000000460.17", 
#"ENSG00000000938.13", 
#"ENSG00000000971.16"))

#cancer_names_assigned <- data.frame(gene_id = formatting_output[["formatted_gene_list"]], gene_name = formatting_output[["gene_names_assigned"]])
#cancer_names_assigned

#print(cancer_names_assigned[cancer_names_assigned$gene_id == "ENSG00000000419.13", "gene_name"])