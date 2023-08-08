setwd("C:/Users/richarej/tcga/tcga_app_07.08.23")

## Loading libraries ##
if (!require("ggcorrplot", quietly = TRUE)) { install.packages("ggcorrplot") } 
if (!require("factoextra", quietly = TRUE)) { install.packages("factoextra") } 
if (!require("ggfortify", quietly = TRUE)) { install.packages("ggfortify") } 
if (!require("umap", quietly = TRUE)) { install.packages("umap") } 
if (!require("plotly", quietly = TRUE)) { install.packages("plotly") } 
if (!require("viridis", quietly = TRUE)) { install.packages("viridis") } 
#devtools::install_github('ropensci/plotly')

# General
library(shiny)
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

## Data ##
source("PCA_function.R")

list_of_cancer_types <- c("Breast Invasive Carcinoma" = "TCGA-BRCA",                                                                             
                          "Thyroid Carcinoma" = "TCGA-THCA",                                                                                     
                          "Uterine Corpus Endometrial Carcinoma" = "TCGA-UCEC",                                                                  
                          "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma" = "TCGA-DLBC",                                                      
                          "Colon Adenocarcinoma" = "TCGA-COAD",                                                                                  
                          "Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma" = "TCGA-CESC",                                      
                          "Bladder Urothelial Carcinoma" = "TCGA-BLCA",                                                                          
                          "Cholangiocarcinoma" = "TCGA-CHOL",                                                                                    
                          "Esophageal Carcinoma" = "TCGA-ESCA",                                                                                  
                          "Adrenocortical Carcinoma" = "TCGA-ACC",                                                                              
                          "Kidney Chromophobe" = "TCGA-KICH",                                                                                    
                          "Head and Neck Squamous Cell Carcinoma" = "TCGA-HNSC",                                                
                          "Liver Hepatocellular Carcinoma" = "TCGA-LIHC",      
                          "Mesothelioma" = "TCGA-MESO",                                                                                         
                          "Acute Myeloid Leukemia" = "TCGA-LAML",                                                                                
                          "Kidney Renal Papillary Cell Carcinoma" = "TCGA-KIRP",                                                                 
                          "Kidney Renal Clear Cell Carcinoma" = "TCGA-KIRC",                                                                     
                          "Glioblastoma Multiforme" = "TCGA-GBM",                                                                               
                          "Brain Lower Grade Glioma" = "TCGA-LGG",                                                                              
                          "Sarcoma" = "TCGA-SARC",                                                                                               
                          "Pheochromocytoma and Paraganglioma" = "TCGA-PCPG",                                                                    
                          "Rectum Adenocarcinoma" = "TCGA-READ",                                                                                 
                          "Pancreatic Adenocarcinoma" = "TCGA-PAAD",                                                                             
                          "Lung Adenocarcinoma" = "TCGA-LUAD",                                                                                   
                          "Prostate Adenocarcinoma" = "TCGA-PRAD",                                                                               
                          "Ovarian Serous Cystadenocarcinoma" = "TCGA-OV",                                                                     
                          "Lung Squamous Cell Carcinoma" = "TCGA-LUSC",                                                                          
                          "Testicular Germ Cell Tumors" = "TCGA-TGCT",                                                                           
                          "Thymoma" = "TCGA-THYM",                                                                                               
                          "Uveal Melanoma" = "TCGA-UVM",                                                                                       
                          "Skin Cutaneous Melanoma" = "TCGA-SKCM",                                                                               
                          "Uterine Carcinosarcoma" = "TCGA-UCS",                                                                                
                          "Stomach Adenocarcinoma" = "TCGA-STAD")


# UI----------------------------------------------------------------------------
ui <- fluidPage(
  
  titlePanel("TCGA Dummy App"),
  
  sidebarLayout(
    sidebarPanel(
      radioButtons("cancer_type_list", "Select Cancer Type", list_of_cancer_types),
      
      textAreaInput("gene_string", "Input your gene names (separated by new lines)", rows = 5), 
      
      actionButton("button_display_plot", "Generate PCA Plots")
    ), # Sidebar panel
    
    mainPanel(
      plotlyOutput("display_plot", width = "100%",
                   height = "1000px")
    ) # Main panel
  ) # Sidebar layout
) # Fluid page


# Server------------------------------------------------------------------------
server <- function(input, output) {
  
  observeEvent(input$button_display_plot, {
    output$display_plot <- renderPlotly({
      pca_function(input$cancer_type_list, input$gene_string)
    }) # Render Plotly
  }) # Observe event
  
} # Server


# Run---------------------------------------------------------------------------
shinyApp(ui = ui, server = server)

#Test genes
#ENSG00000000003.15
#ENSG00000000005.6
