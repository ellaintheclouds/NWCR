# Set up------------------------------------------------------------------------
setwd("D:/app versions/tcga_app_usb")

## Loading libraries ##
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")

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

# UI----------------------------------------------------------------------------
ui <- fluidPage(
  
  titlePanel("TCGA Dummy App"),
  
  sidebarLayout(
    sidebarPanel(
      checkboxGroupInput("cancer_type_list", "Select Cancer Type", list_of_cancer_types),
      br(), 
      
      print(HTML("<strong>Input gene names</strong> <br/> use Ensembl format and separate by new lines")), 
      textAreaInput("gene_string", "", rows = 5), 
      
      actionButton("button_pca_analysis", "Principal Component Analysis"), 
      
    ), # Sidebar panel
    
    mainPanel(
      uiOutput("cancer_type_menu"), 
      uiOutput("gene_menu"), 
      uiOutput("pcx_menu"), 
      uiOutput("pcy_menu"), 
      
      plotlyOutput("display_plot", width = "100%",
                   height = "1000px")
    ) # Main panel
  ) # Sidebar layout
) # Fluid page


# Server------------------------------------------------------------------------
server <- function(input, output) {
  
  observeEvent(input$button_pca_analysis, {
    
    out_plots <- pca_function(input$cancer_type_list, input$gene_string)
    
    # Display plot (interactive)-----
    # Setting default variables
    first_cancer_type <- input$cancer_type_list[[1]]
    
    gene_list <- strsplit(input$gene_string, split = "\n")[[1]] # Split the user input into a list
    first_gene <-  gene_list[1]
    pcx <- 1
    pcy <- 2
    
    # UI variable menus
    output$cancer_type_menu <- renderUI({
      selectInput("display_cancer_type", "Cancer type", input$cancer_type_list)
    })
    
    output$gene_menu <- renderUI({
      selectInput("display_gene", "Gene", gene_list)
    })
    
    output$pcx_menu <- renderUI({
      selectInput("display_pcx", "PC (x-axis)", c(1:10))
    })
    
    output$pcy_menu <- renderUI({
      selectInput("display_pcy", "PC (y-axis)", c(1:10))
    })
    
    # Plot
    output$display_plot <- renderPlotly({
      out_plots[[first_cancer_type]][[first_gene]][[paste0(pcx, "_", pcy)]]
      
      #--------------------------------------------------------------------------------------------------------------------------------------
      #      
      #    # Download display plot
      #      output$download_display = downloadHandler(
      #        filename = function() {"plots.pdf"},
      #        content = function(file) {
      #          pdf(file, onefile = TRUE, width = 15, height = 9)
      #          replayPlot(out_plots[[first_cancer_type]][[first_gene]][[paste0(pcx, "_", pcy)]])
      #          dev.off()
      #        } # Function
      #      ) # Download handler
      #      
      #    # Download all plots-----
      #    output$download_all = downloadHandler(
      #      filename = function() {"plots.pdf"},
      #      content = function(file) {
      #        pdf(file, onefile = TRUE, width = 15, height = 9)
      #        replayPlot(out_plots[[first_cancer_type]][[first_gene]][[paste0(pcx, "_", pcy)]])
      #        dev.off()
      #      } # Function
      #    ) # Download handler
      #      
      #--------------------------------------------------------------------------------------------------------------------------------------
      
    }) # Render Plotly
  }) # Observe event
  
} # Server


# Run---------------------------------------------------------------------------
shinyApp(ui = ui, server = server)

#Test genes
#ENSG00000000003.15
#ENSG00000000005.6
