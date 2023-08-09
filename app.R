# Set up------------------------------------------------------------------------
setwd("D:/app versions/tcga_app_usb")

## Loading libraries ##
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("TCGAbiolinks")

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("edgeR")

if (!require("ggcorrplot", quietly = TRUE)) { install.packages("ggcorrplot") } 
if (!require("factoextra", quietly = TRUE)) { install.packages("factoextra") } 
if (!require("ggfortify", quietly = TRUE)) { install.packages("ggfortify") } 
if (!require("umap", quietly = TRUE)) { install.packages("umap") } 
if (!require("plotly", quietly = TRUE)) { install.packages("plotly") } 
if (!require("viridis", quietly = TRUE)) { install.packages("viridis") } 
if (!require("zip", quietly = TRUE)) { install.packages("zip") } 

# General
library(shiny)
library(stringr) # Makes working with strings simpler.
library(viridis) # Colour scheme.
library(zip)

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

# Data----------
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
  
  # Sidebar----------
  sidebarLayout(
    sidebarPanel(
      # Cancer type input
      checkboxGroupInput("cancer_type_list", "Select Cancer Type", list_of_cancer_types),
      br(), 
      
      # Gene input
      print(HTML("<strong>Input gene names</strong> <br/> use Ensembl format and separate by new lines")), 
      textAreaInput("gene_string", "", rows = 5), 
      
      # Button that triggers PCA computation
      actionButton("button_pca_analysis", "Principal Component Analysis"), 
      
    ), # Sidebar panel
    
    # Main panel----------
    mainPanel(
      # Drop-down menus that reactively alter display PCA plot
      uiOutput("cancer_type_menu"), 
      uiOutput("gene_menu"), 
      uiOutput("pcx_menu"), 
      uiOutput("pcy_menu"), 

      # Display PCA plot
      plotlyOutput("display_plot", width = "100%",
                   height = "1000px"), 
      
      # Download buttons
      uiOutput("download_pdf_button"), 
      uiOutput("download_zip_button")

      
    ) # Main panel
  ) # Sidebar layout
) # Fluid page


# Server------------------------------------------------------------------------
server <- function(input, output) {
  
  observeEvent(input$button_pca_analysis, {
    
    out_plots <- pca_function(input$cancer_type_list, input$gene_string)
    
    # Display plot (interactive)-----
    # UI variable menus (with the first variable as default
    output$cancer_type_menu <- renderUI({
      selectInput("display_cancer_type", "Cancer type", input$cancer_type_list, selected = input$cancer_type_list[1])
    })
    
    output$gene_menu <- renderUI({
      selectInput("display_gene", "Gene", strsplit(input$gene_string, split = "\n")[[1]], selected = strsplit(input$gene_string, split = "\n")[[1]][1])
    })
    
    output$pcx_menu <- renderUI({
      selectInput("display_pcx", "PC (x-axis)", c(1:10), selected = 1)
    })
    
    output$pcy_menu <- renderUI({
      selectInput("display_pcy", "PC (y-axis)", c(1:10), selected = 2)
    })
    
    # Selected data (reactive)
    display_cancer_type_data <- reactive({input$display_cancer_type})
    display_gene_data <- reactive({input$display_gene})
    display_pcx_data <- reactive({input$display_pcx})
    display_pcy_data <- reactive({input$display_pcy})
    
    # Plot
    output$display_plot <- renderPlotly({
      out_plots[[input$display_cancer_type]][[input$display_gene]][[paste0(input$display_pcx, "_", input$display_pcy)]]
    }) # Render Plotly
    
    # Download all plots-------------------------------------------------------------------------------------------
    # Reactive download buttons
    output$download_zip_button <- renderUI(downloadButton("download_zip", "Download All Plot Images"))
    output$download_pdf_button <- renderUI(downloadButton("download_pdf", "Download All Plots as a PDF"))
    
    # Create .zip with individual .png plots
    output$download_zip <-downloadHandler(
      filename = function(){"out.zip"},
      content = function(file){
        if(TRUE){
        # create png folder
        dir.create("out_download", showWarnings = FALSE)
        # loop through cancer type
        for(current_cancer_type in names(out_plots)){
          dir.create(file.path(paste0("out_download/", current_cancer_type)), showWarnings = FALSE)
          # loop through genes
          for(current_gene in names(out_plots[[current_cancer_type]])){
            dir.create(file.path(paste0("out_download/", current_cancer_type, "/", current_gene)), showWarnings = FALSE)
            # loop through pca combinations
            for(pca_combination in names(out_plots[[current_cancer_type]][[current_gene]])){
              # save each graph
              out_plot_fp <- paste0("out_download/", current_cancer_type, "/", current_gene, "/", pca_combination, ".png")
              ggsave(out_plot_fp, out_plots[[current_cancer_type]][[current_gene]][[pca_combination]])
            }
          }        
        }
        all_files <- list.files("out_download", full.names = TRUE)
        }

        zip::zip(zipfile = file, files=all_files)
      }, # Content function
      
      contentType = "application/zip"
      ) # Download handler
    
    
    # Create singular .pdf containing all plots
    output$download_pdf = downloadHandler(
      filename = function() {"plots.pdf"}, 
      content = function(file) {
        pdf(file, onefile = TRUE, width = 15, height = 10)
        for(plot_printout in out_plots){ print(plot_printout)}
        dev.off()
      })                      
    
  }) # Observe event
} # Server


# Run---------------------------------------------------------------------------
shinyApp(ui = ui, server = server)

#Test genes ENSG00000000003.15ENSG00000000005.6