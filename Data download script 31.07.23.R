### Set up ###
  setwd("C:/Users/richarej/tcga")
  
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("TCGAbiolinks")
  
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("summarizedExperiment")
  
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("edgeR")
  
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(edgeR)


### Providing a list of the cancer types
  cancer_codes_passed <- list()
  
  cancer_codes <- c("TCGA-BRCA", "TCGA-THCA", "TCGA-UCEC", "TCGA-DLBC", "TCGA-COAD", 
                    "TCGA-CESC", "TCGA-BLCA", "TCGA-CHOL", "TCGA-ESCA", "TCGA-ACC", 
                    "TCGA-KICH", "TCGA-HNSC", "TCGA-LIHC", "TCGA-MESO", "TCGA-LAML", 
                    "TCGA-KIRP", "TCGA-KIRC", "TCGA-GBM", "TCGA-LGG", "TCGA-SARC", 
                    "TCGA-PCPG", "TCGA-READ", "TCGA-PAAD", "TCGA-LUAD", "TCGA-PRAD", 
                    "TCGA-OV", "TCGA-LUSC", "TCGA-TGCT", "TCGA-THYM", "TCGA-UVM", 
                    "TCGA-SKCM", "TCGA-UCS", "TCGA-STAD")
  
  # These projects possibly have errors. We need to come back to these but will exclude them for now.
  ignore_project <- c("TCGA-UCEC", "TCGA-BLCA", "TCGA-PAAD", "TCGA-KIRC", "TCGA-KIRP", "TCGA-MESO")
  
  retry_projects <- list()
  less_sample_projects <- list()
  num_sample <- list()

### Query ###
for(x in cancer_codes){
  project <- x
  # If the current project ID is not in ignore_project or cancer_codes_passed (from already succeeding), run query:
  if( !project %in% ignore_project & ! project %in% names(cancer_codes_passed)){
    query_RNASeq <- tryCatch({GDCquery(project = project , # Try  
                                       data.category = "Transcriptome Profiling", 
                                       data.type = "Gene Expression Quantification", 
                                       workflow.type = "STAR - Counts", 
                                       experimental.strategy = "RNA-Seq" )}, 
                             error =  function(cond){
                                      message(cond) # Provide an error message
                                      return(NULL)},
                             warning = function(cond){
                                      message(cond) # Provide a warning message
                                      return(NULL)})
    RNAseq_samples <- tryCatch({ getResults(query_RNASeq)}, 
                             error = function(cond){return(NULL)})   
    
    if(!is.null(query_RNASeq)){
      if(!is.null(RNAseq_samples)){
        cancer_codes_passed[[x]] <- query_RNASeq
        num_sample[[x]] <- nrow(RNAseq_samples)
        
      } else {
        less_sample_projects[[project]] = NULL
      }
    }
  }
}


### Download ###
for (x in cancer_codes_passed){
  GDCdownload(x)  
}
