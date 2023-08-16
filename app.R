# Set up-------------------------------------------------------------------------
# Directory----------
setwd("D:/app versions/tcga_app_usb")

# Libraries----------
if (!require("bslib", quietly = TRUE)) { install.packages("bslib") } 
if (!require("shinythemes", quietly = TRUE)) { install.packages("shinythemes") } 
if (!require("ggcorrplot", quietly = TRUE)) { install.packages("ggcorrplot") } 
if (!require("factoextra", quietly = TRUE)) { install.packages("factoextra") } 
if (!require("ggfortify", quietly = TRUE)) { install.packages("ggfortify") } 
if (!require("umap", quietly = TRUE)) { install.packages("umap") } 
if (!require("plotly", quietly = TRUE)) { install.packages("plotly") } 
if (!require("viridis", quietly = TRUE)) { install.packages("viridis") } 
if (!require("zip", quietly = TRUE)) { install.packages("zip") } 
if (!require("edgeR", quietly = TRUE)) { BiocManager::install("edgeR") } 
if (!require("TCGAbiolinks", quietly = TRUE)) { BiocManager::install("TCGAbiolinks") } 

# General
library(shiny)
library(stringr) # Makes working with strings simpler.
library(viridis) # Colour scheme.
library(zip) # Allows creation of zip files

# UI
library(shinythemes)
library(forcats)
library(plyr)
library(shinyWidgets)

# Biological data
library(TCGAbiolinks) # For integrative analysis of GDC data.
library(ggplot2)
library(ggcorrplot) # Allows visulaisation of a correlation matrix using ggplot2.
library(SummarizedExperiment) # For storing processed data from high-throughput sequencing assays.
library(edgeR) # Differential expression analysis of RNA-seq expression profiles.
library(plotly) # Creates interactive, publication-suitable graphs.

# PCA
library(factoextra) # For principal component analysis
library(ggfortify) # Allows plotting of PCA and survival analysis.
library(umap) # Algorithm for dimensional reduction

# Data----------
source("PCA_function.R")
source("input_name_function.R")

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

list_of_cancer_types_rev <- vector()
for(current_name in names(list_of_cancer_types)){
  current_code <- list_of_cancer_types[current_name]
  list_of_cancer_types_rev[current_code] = current_name
}


# UI----------------------------------------------------------------------------
ui <- fluidPage(theme = shinytheme("sandstone"), 
                
                # Setting the font for the entire page
                tags$head(
                  tags$head(tags$style(HTML('* {font-family: "Calibri"};')))
                ),
                
                navbarPage("TCGA Survival and Principal Component Analysis", 
                           
                           # Analysis-------------------------------------------
                           tabPanel("Analysis", 
                                    
                                    sidebarLayout(
                                      sidebarPanel( width = 3,
                                                    # Gene input
                                                    print(HTML("<strong>Input gene names</strong> <br/> format as gene name or Ensembl ID <br/> (separate by new lines) ")), 
                                                    textAreaInput("gene_string", "", rows = 5), 
                                                    br(), 
                                                    
                                                    # Cancer type input
                                                    selectInput("cancer_type_list", "Select Cancer Type", list_of_cancer_types, multiple = TRUE),
                                                    
                                                    # Compute button
                                                    actionButton("button_compute", "Analyse"),
                                                    br(), 
                                                    
                                                    # Reactive error message
                                                    textOutput("invalid_gene_message_display")
                                                    
                                      ), # Sidebar panel 
                                      
                                      mainPanel(
                                        tabsetPanel(type = "tabs", 
                                                    
                                                    # Principal component analysis-----------------------------------------
                                                    tabPanel("Principal Component Analysis",
                                                             
                                                             fluidRow(
                                                               # Drop-down menus----------
                                                               column(width = 3, uiOutput("cancer_type_menu")), 
                                                               column(width = 3, uiOutput("gene_menu")), 
                                                               column(width = 3, uiOutput("pcx_menu")), 
                                                               column(width = 3, uiOutput("pcy_menu"))
                                                             ), # Fluid row
                                                             
                                                             fluidRow(
                                                               # Display PCA plot----------
                                                                      plotlyOutput("display_pca_plot", width = "100%", height = "500px")
                                                             ), # Fluid row
                                                               
                                                            fluidRow(
                                                               # Display scree plot----------
                                                                      plotlyOutput("display_scree_plot", width = "100%", height = "300px")
                                                             ), # Fluid row
                                                             
                                                             fluidRow(
                                                               # Download buttons----------
                                                               column(width = 6, uiOutput("download_pdf_button")), 
                                                               column(width = 6, uiOutput("download_zip_button"))
                                                             ) # Fluid row
                                                    ), # Tab panel
                                                    
                                                    # Survival analysis------------------------------
                                                    tabPanel("Survival Analysis")
                                        ) # Tabset panel
                                      ) # Main panel
                                    ) # Sidebar layout
                           ), # Tab panel

                           # About----------------------------------------------
                           tabPanel("About", 
                                    fluidRow(
                                      
                                      column(width = 4, wellPanel(
                                        HTML("Created by A. Kennedy, E. Richardson and B. Shih. <br/> <br/> 
                      This tool was created with funding from North West Cancer Research.")
                                      ) # Well panel
                                      ), # Column
                                      
                                      column(width = 8, wellPanel(
                                        HTML('<iframe width="560" height="315" src="https://www.youtube.com/embed/Ka2pWqXS1WA" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen></iframe>')
                                      ) # Well panel
                                      ) # Column
                                    ) # Fluid row
                           ), # Tab panel
                           
                           # Contact us-------------------------------------------------------
                           tabPanel("Contact us"), 
                           
                           # Source code------------------------------------------------------
                           tabPanel("Source code", 
                                    HTML('<script src="https://emgithub.com/embed-v2.js?target=https%3A%2F%2Fgithub.com%2Fellaintheclouds%2FNWCR%2Fblob%2Fmain%2FPCA_function.R&style=default&type=code&showBorder=on&showLineNumbers=on&showFileMeta=on&showFullPath=on&showCopy=on"></script>')
                           )
                ) # Navbar page
) # Fluid page


# Server------------------------------------------------------------------------
server <- function(input, output, session){
  observeEvent(input$button_compute, {
    
    # Formatting function----------
    gene_list <- strsplit(input$gene_string, split = "\n")[[1]] # Split the user input into a list
    formatting_output <- input_format(gene_list)
    print(head(formatting_output[["formatted_gene_list"]]))
    formatted_gene_list <- formatting_output[["formatted_gene_list"]]$gene_id
    names(formatted_gene_list) <- formatting_output[["formatted_gene_list"]]$merged_name
    print(formatted_gene_list)
    
    invalid_gene_message <- reactive({
        if (length(formatting_output[["rejected_gene_list"]]) > 0){
        paste0("These inputs are not valid gene names: ", formatting_output[["rejected_gene_list"]], 
               ". Make sure that the format, letter case, and line separation is correct.")
        } # If
      else {""}
    }) # Reactive
    
    output$invalid_gene_message_display <- renderText({invalid_gene_message()})
    
    # PCA function----------
    output_data <- pca_function(input$cancer_type_list, formatted_gene_list)
    
    # Re-assigning names to cancer types
    input_cancer_type <- input$cancer_type_list
    names(input_cancer_type) <- list_of_cancer_types_rev[input_cancer_type]
    
    # Display plot (interactive)----------
    # UI variable menus (with the first variable as default
    print("DONE, ABOUT TO MAKE BUTTON")
    
    output$cancer_type_menu <- renderUI({
      selectInput("display_cancer_type", "Cancer type", input_cancer_type, selected = input_cancer_type[1])
    })
    
    output$gene_menu <- renderUI({
      selectInput("display_gene", "Gene", formatted_gene_list, selected = formatted_gene_list[1])
    })
    
    output$pcx_menu <- renderUI({
      selectInput("display_pcx", "PC (x-axis)", c(1:9), selected = 1)
    })
    
    # Making sure that only valid y inputs are displayed (relative to x)
    observeEvent(input$display_pcx, {
      
      if (input$display_pcx == 1){pcy_choices <- c(2:10)}
      else if (input$display_pcx == 2){pcy_choices <- c(3:10)}
      else if (input$display_pcx == 3){pcy_choices <- c(4:10)}
      else if (input$display_pcx == 4){pcy_choices <- c(5:10)}
      else if (input$display_pcx == 5){pcy_choices <- c(6:10)}
      else if (input$display_pcx == 6){pcy_choices <- c(7:10)}
      else if (input$display_pcx == 7){pcy_choices <- c(8:10)}
      else if (input$display_pcx == 8){pcy_choices <- c(9:10)}
      else if (input$display_pcx == 9){pcy_choices <- 10}
      
      output$pcy_menu <- renderUI({selectInput("display_pcy", "PC (y-axis)", pcy_choices, selected = pcy_choices[1])})
    })
    
    # Retrieving PCA plot----------
    output$display_pca_plot <- renderPlotly({
      validate(need(input$display_cancer_type, input$display_gene, message = FALSE)) # Validate needs
      output_data[["pca plots"]][[input$display_cancer_type]][[input$display_gene]][[paste0(input$display_pcx, "_", input$display_pcy)]]
    }) # Render Plotly
    
    # Plotting scree plot----------
    output$display_scree_plot <- renderPlotly({
      
      # Validate needs
      validate(need(input$display_cancer_type, message = FALSE))
      
      # Get data
      screeplot_df <- output_data[["contribution percentile dataframes"]][[input$display_cancer_type]]
      current_gene_name <- formatting_output[["formatted_gene_list"]][formatting_output[["formatted_gene_list"]]$gene_id == input$display_gene, "gene_name"]
      
      # Plotting for expressed genes
      if(input$display_gene %in% colnames(screeplot_df)){
        current_scree_data <- screeplot_df[,c("PC", "contribution", input$display_gene)]
        colnames(current_scree_data) <- c("PC", "contribution", "current_gene_column") 
        
        percent_label <- paste0(current_gene_name, "\nContribution\n(%)")
        
        scree_plot <- ggplot(data = current_scree_data, mapping = aes(x = PC, y = contribution, fill = current_gene_column)) +
          geom_bar(stat = "identity") +
          labs(fill = percent_label) + xlab("Principal Component") + ylab("Variance Explained by PC (%)") + 
          scale_x_discrete(limit = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")) +
          theme_minimal() + 
          scale_fill_viridis(option = "magma")
        
        # Plotting for genes with almost no expression
      } else {
        current_scree_data <- screeplot_df[,c("PC", "contribution")]
        colnames(current_scree_data) <- c("PC", "contribution") 
        
        scree_plot <- ggplot(data = current_scree_data, mapping = aes(x = PC, y = contribution)) +
          geom_bar(stat = "identity") + xlab("Principal Component") + ylab("Variance Explained by PC (%)") + 
          scale_x_discrete(limit = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")) +
          theme_minimal() + 
          scale_fill_viridis(discrete = TRUE, option = "magma")
      } # Else   
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
          for(current_cancer_type in names(output_data[["pca plots"]])){
            dir.create(file.path(paste0("out_download/", current_cancer_type)), showWarnings = FALSE)
            # loop through genes
            for(current_gene in names(output_data[["pca plots"]][[current_cancer_type]])){
              dir.create(file.path(paste0("out_download/", current_cancer_type, "/", current_gene)), showWarnings = FALSE)
              # loop through pca combinations
              for(pca_combination in names(output_data[["pca plots"]][[current_cancer_type]][[current_gene]])){
                # save each graph
                out_plot_fp <- paste0("out_download/", current_cancer_type, "/", current_gene, "/", pca_combination, ".png")
                ggsave(out_plot_fp,
                       plot = output_data[["pca plots"]][[current_cancer_type]][[current_gene]][[pca_combination]], 
                       width = 8, height = 6
                ) # GGsave
              } # For
            } # For        
          } # For
          all_files <- list.files("out_download", full.names = TRUE)
        } # If
        
        zip::zip(zipfile = file, files = all_files)
      }, # Content function
      
      contentType = "application/zip"
    ) # Download handler
    
    # Create singular .pdf containing all plots
    output$download_pdf = downloadHandler(
      filename = function() {"plots.pdf"}, 
      content = function(file) {
        pdf(file, onefile = TRUE, width = 8, height = 6)
        for(plot_printout in output_data[["pca plots"]]){ print(plot_printout)}
        dev.off()
      } # Content function
    ) # Download handler                      
  }) # Observe event
} # Server


# Run---------------------------------------------------------------------------
shinyApp(ui = ui, server = server)
