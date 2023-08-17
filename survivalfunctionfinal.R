getwd()
list.files()
#### load libraries 
if (!require("ggcorrplot", quietly = TRUE)) { install.packages("ggcorrplot") } 
if (!require("viridis", quietly = TRUE)) { install.packages("viridis") }
if (!require("zip", quietly = TRUE)) { install.packages("zip") } 
library(ggcorrplot)
library(stringr) 
library(viridis)
library(zip)
library(shiny) 
library(TCGAbiolinks)
library(edgeR)
library(SummarizedExperiment)
library(here)

####
#source in function
source("modSurvKM.r")


#list cancer types
list_of_cancer_types <- c("Acute Myeloid Leukemia" = "TCGA-LAML",                                                                                
                          "Adrenocortical Carcinoma" = "TCGA-ACC",                                                                              
                          "Brain Lower Grade Glioma" = "TCGA-LGG",                                                                              
                          "Breast Invasive Carcinoma" = "TCGA-BRCA",                                                                             
                          "Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma" = "TCGA-CESC",                                      
                          "Cholangiocarcinoma" = "TCGA-CHOL",                                                                                    
                          "Colon Adenocarcinoma" = "TCGA-COAD",                                                                                  
                          "Esophageal Carcinoma" = "TCGA-ESCA",                                                                                  
                          "Glioblastoma Multiforme" = "TCGA-GBM",                                                                               
                          "Head and Neck Squamous Cell Carcinoma" = "TCGA-HNSC",                                                
                          "Kidney Chromophobe" = "TCGA-KICH",                                                                                    
                          "Liver Hepatocellular Carcinoma" = "TCGA-LIHC",      
                          "Lung Adenocarcinoma" = "TCGA-LUAD",                                                                                   
                          "Lung Squamous Cell Carcinoma" = "TCGA-LUSC",                                                                          
                          "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma" = "TCGA-DLBC",                                                      
                          "Ovarian Serous Cystadenocarcinoma" = "TCGA-OV",                                                                     
                          "Pheochromocytoma and Paraganglioma" = "TCGA-PCPG",                                                                    
                          "Prostate Adenocarcinoma" = "TCGA-PRAD",                                                                               
                          "Rectum Adenocarcinoma" = "TCGA-READ",                                                                                 
                          "Sarcoma" = "TCGA-SARC",                                                                                               
                          "Skin Cutaneous Melanoma" = "TCGA-SKCM",                                                                               
                          "Stomach Adenocarcinoma" = "TCGA-STAD", 
                          "Testicular Germ Cell Tumors" = "TCGA-TGCT",                                                                           
                          "Thymoma" = "TCGA-THYM",                                                                                               
                          "Thyroid Carcinoma" = "TCGA-THCA",                                                                                     
                          "Uterine Carcinosarcoma" = "TCGA-UCS",                                                                                
                          "Uveal Melanoma" = "TCGA-UVM")

## create vector of cancer types 


list_of_cancer_types_fun <- vector()
for(current_name in names(list_of_cancer_types)){
  current_cancer_type <- list_of_cancer_types[current_name]
  list_of_cancer_types_fun[current_cancer_type] = current_name
}

##create function START
survival_function <- function(cancer_type_list, gene_list){ 
  
  
  ###string split the list of genes that user input
 

  
  
  ###store the plots as a list 
  output_survival_plot <- list()

  
  ### loop through the cancer type data and give names 
  for(current_cancer_type in cancer_type_list){
    
    ### reassign names 
    #current_cancer_type_named <- current_cancer_type 
    names(current_cancer_type) <- list_of_cancer_types_fun[current_cancer_type]
    
    current_file_name <- paste0(current_cancer_type, ".rds")
    print(current_file_name)
    
   output_survival_plot[[current_cancer_type]] <- list()
    
    ### 
    # read in RDS files and run function 
    
    ### read in RDS file of cancer type that will be used in loop
    cancer_data <- readRDS(paste0("data/", current_file_name))    
    # query to get clinical data for the current cancer project
    data_clinic <- GDCquery_clinic(current_cancer_type, "clinical")
    count_mx <- assay(cancer_data, "unstranded")
   gene_annotation <- rowData(cancer_data)
    dgelist <- DGEList(counts=count_mx, genes=gene_annotation)
    dgelist <- calcNormFactors(dgelist)
    normcounts <- cpm(dgelist, normalized.lib.sizes=T)
    tmp_reorganised_log2 <-log2(normcounts +1)
    print(gene_list)
  
    tabSurvKM <- TCGAanalyze_SurvivalKM(
      data_clinic,
      tmp_reorganised_log2,
      Genelist= gene_list,
      Survresult= T,
      p.cut =1,
      ThreshTop = 0.67,
      ThreshDown = 0.33
    )
    output_survival_plot[[current_cancer_type]] <-tabSurvKM
    
  }
  
  return(output_survival_plot)
  
}

