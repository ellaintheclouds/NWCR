setwd("C:/Users/richarej/tcga//app versions/tcga_app_1_07.08.23")

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

list_of_cancer_types <- c("Acute Myeloid Leukemia" = "TCGA-LAML",                                                                                
                          "Adrenocortical Carcinoma" = "TCGA-ACC",                                                                              
                          "Bladder Urothelial Carcinoma" = "TCGA-BLCA",                                                                          
                          "Brain Lower Grade Glioma" = "TCGA-LGG",                                                                              
                          "Breast Invasive Carcinoma" = "TCGA-BRCA",                                                                             
                          "Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma" = "TCGA-CESC",                                      
                          "Cholangiocarcinoma" = "TCGA-CHOL",                                                                                    
                          "Colon Adenocarcinoma" = "TCGA-COAD",                                                                                  
                          "Esophageal Carcinoma" = "TCGA-ESCA",                                                                                  
                          "Glioblastoma Multiforme" = "TCGA-GBM",                                                                               
                          "Head and Neck Squamous Cell Carcinoma" = "TCGA-HNSC",                                                
                          "Kidney Chromophobe" = "TCGA-KICH",                                                                                    
                          "Kidney Renal Clear Cell Carcinoma" = "TCGA-KIRC",                                                                     
                          "Kidney Renal Papillary Cell Carcinoma" = "TCGA-KIRP",                                                                 
                          "Liver Hepatocellular Carcinoma" = "TCGA-LIHC",      
                          "Lung Adenocarcinoma" = "TCGA-LUAD",                                                                                   
                          "Lung Squamous Cell Carcinoma" = "TCGA-LUSC",                                                                          
                          "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma" = "TCGA-DLBC",                                                      
                          "Mesothelioma" = "TCGA-MESO",                                                                                         
                          "Ovarian Serous Cystadenocarcinoma" = "TCGA-OV",                                                                     
                          "Pancreatic Adenocarcinoma" = "TCGA-PAAD",                                                                             
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
                          "Uterine Corpus Endometrial Carcinoma" = "TCGA-UCEC",                                                                  
                          "Uveal Melanoma" = "TCGA-UVM")


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
    out_plots <- pca_function(input$cancer_type_list, input$gene_string)
    first_cancer_type <- input$cancer_type_list
    first_gene <-  strsplit(input$gene_string, split = "\n")[[1]][1]
    pcx <- 1
    pcy <- 2
    
    output$display_plot <- renderPlotly({
      out_plots[[first_cancer_type]][[first_gene]][[paste0(pcx, "_", pcy)]]
    }) # Render Plotly
  }) # Observe event
  
} # Server


# Run---------------------------------------------------------------------------
shinyApp(ui = ui, server = server)

#Test genes
#ENSG00000000003.15
#ENSG00000000005.6
