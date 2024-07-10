library(Seurat)
library(shiny)
library(DT)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(pals)
library(dplyr)
library(bslib)
library(pheatmap)
library(DESeq2)

debug <- T

ui <- navbarPage(
  title = "HS-OmicsDB: Hidradenitis Suppurativa Omics Database",
  tabPanel(
    title = "Home",
    br(),
    hr(),
    p(
      style = "text-align: justify; font-size = 25px",
      "Here we present HS 'omics database (HS-OmicsDB)',
          an initiative to integrate the growing body of 'omics data'
          collected from studies of HS. We hope this initiative provides
          a tool for the HS research community to enhance investigations
          into the etiology of and specific functions of novel disease driving genes.
          Our database compiles publicly available high-throughput bulk-RNA and single
          cell RNA sequencing studies with publicly available fastq files from studies
          on patients with HS deposited on the National institute of Health gene expression omnibus (GEO)."
    ),
    p("Users can explore genes of interest to determine tissue-specific and cell-specific gene expression: "),
    p("Current data from bulk-sequencing studies includes: "),
    tags$ul(
      tags$li("73 lesional (L)"),
      tags$li("28 perilesional (PL)"),
      tags$li("30 nonlesional (NL)"),
      tags$li("19 healthy control (HC)")
    ),
    p("Current single cell sequencing data includes 11 samples from lesional HS skin"),
    p("With the wave of thriving HS research, this initiative could serve as a valuable
          resource to increase access to valuable sequencing data. HS-OmicsDB will prove to be
          a valuable validation and hypothesis generating tool in the study of HS. Future versions
          will incorporate the growing number of ongoing HS studies employing 'omics technologies'
          (e.g. proteomics, lipidomics, mass cytometry, spatial sequencing etc.)
          to further enhance research in the HS community. We hope to build this resource to
          provide a 'Rosetta stone' of the complex cutaneous dysregulation that occurs
          in HS by integrating data from multiple 'omics technologies',
          that is easy to be use and flexible in its application.
          Such a tool will help build an integrative view of the dynamic process
          that occurs throughout the progression of HS. For interest in getting
          involved with the HSOmicsDB Initiative please contact Qing-Sheng Mi (qmi1@hfhs.org)."),
    p("Bulk RNAseq Datasets:"),
    p("GSE151243, GSE154773, GSE155176, GSE189266, GSE213761"),
    p("Single-cell RNAseq Datasets:"),
    p("GSE154775, GSE175990"),
    hr()
  ),
  tabPanel(
    title = "Visualize gene expression",
    fluidPage(
      sidebarLayout(
        sidebarPanel(
          width = 3,
          h3("Data"),
          textInput("gene", "Genes", "IL17A"),
          textOutput("warning_message"),
          actionButton("submit", "Submit"),
          h3("Violin plot"),
          numericInput("pt_size", "Point size", value = 0.1, min = 0, max = 1, step = 0.1)
        ),
        mainPanel(
          plotOutput("dr_plot", width = "100%"),
          fluidRow(
            splitLayout(cellWidths = c("50%", "50%"), plotOutput("feature"), plotOutput("vln"))
          ),
          fluidRow(
            splitLayout(cellWidths = c("50%", "50%"), downloadLink("download_feature", "Download PDF"), downloadLink("download_vln", "Download PDF"))
          ),
          plotOutput("dot"),
          downloadLink("download_dot", "Download PDF"),
          # plotOutput("vln"),
          # downloadLink("download_vln", "Download PDF"),
          plotOutput("Correlation")
        )
      )
    )
  ),
  tabPanel(
    title = "Cluster info",
    tableOutput("ncell_table")
  ),
  tabPanel(
    title = "Cluster markers",
    ## p("Use the navigation tools to search, filter, and order the table of gene markers."),
    DT::dataTableOutput("markers")
  ),
  h6("This portal is developed and maintained by the Center for Bioinformatics, Department of Public Health Sciences and the Center for Cutaneous Biology and Immunology Research at Henry Ford Health, Detroit, Michigan.")
)



server <- function(input, output) {
  # reactive will dynamically reload data when RDS is updated
  # no need since will be refreshed each time.

  message("CWD: ", getwd())

  ## TEMP DATA DIR ##
  # removed reactive so only loaded once
  ## Todo: why does it run 5 times?
  data <- reactiveVal(NULL)

  observe({
    if (is.null(data())) {
      data_loading_time <- system.time({
        seurat <- readRDS("../../data/hsomicsdb/seurat/HS_merged_lite.RDS")
        cp <- readRDS("../../data/hsomicsdb/Combat_Adjusted_CPM.rds")
        mk <- read.delim("../../data/hsomicsdb/markers/hs.markers.top30_genes.csv",
          sep = ",", stringsAsFactors = FALSE
        )
        vsd <- readRDS("../../data/hsomicsdb/Updated_VSD.rds")
      })
      data(list(seurat = seurat, cp = cp, mk = mk, vsd = vsd))
      if (debug) {
        message(paste("Loading Data Time: ", data_loading_time, "\n"))
      }
    }
  })

  # check if gene exists in seurat object
  valid_input <- function(gene) {
    if (gene %in% rownames(data()$seurat)) {
      output$warning_message <- renderText({
        ""
      })
      return(TRUE)
    } else {
      output$warning_message <- renderText({
        paste("Gene ", gene, " not found in the dataset")
      })
      return(FALSE)
    }
  }

  generate_plots <- function(gene) {
    # make sure datasets have been loaded
    req(data())
    seurat <- data()$seurat
    cp <- data()$cp
    vsd <- data()$vsd
    mk <- data()$mk

    # UMAP plot
    umap_loading_time <- system.time({
      output$dr_plot <- renderPlot({
        Seurat::DimPlot(seurat,
          reduction = "umap",
          pt.size = 1.1, label = TRUE, cols = as.character(polychrome(34)),
          raster = FALSE, label.size = 6
        ) + NoLegend()
      })
    })

    if (debug) {
      message(paste("UMAP Plot Time: ", umap_loading_time, "\n"))
    }

    # Feature plot
    feature_plot_time <- system.time({
      feature_plot <- Seurat::FeaturePlot(seurat,
        features = input$gene,
        reduction = "umap",
        pt.size = 1, raster = FALSE, cols = c(viridis(5))
      )

      output$feature <- renderPlot({
        feature_plot
      })
    })

    if (debug) {
      message(paste("Feature Plot Time: ", feature_plot_time, "\n"))
    }

    # Feature plot download
    output$download_feature <- downloadHandler(
      filename = "feature.pdf",
      content = function(file) {
        ggsave(filename = file, plot = feature_plot(), width = 5 * length(input$gene), height = 4)
      }
    )

    # Violin plot
    vln_plot_time <- system.time({
      vln_plot <- Seurat::VlnPlot(seurat,
        features = input$gene,
        pt.size = input$point_size
      ) + NoLegend() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

      output$vln <- renderPlot({
        vln_plot
      })
    })

    if (debug) {
      message(paste("Violin Plot Time: ", vln_plot_time, "\n"))
    }

    # Violin plot download
    output$download_vln <- downloadHandler(
      filename = "violin.pdf",
      content = function(file) {
        ggsave(filename = file, plot = vln_plot(), width = 7 * length(input$gene), height = 5)
      }
    )

    # Dot plot
    dot_plot_time <- system.time({
      dot_plot <- ggplot(cp, aes(x = factor(tissue_type, levels = c("Lesional", "Perilesional", "Non-lesional", "Control")), y = eval(as.name(input$gene)), fill = tissue_type)) +
        geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5, stroke = NA) +
        scale_fill_manual(values = c("Lesional" = "firebrick4", "Perilesional" = "darkgoldenrod3", "Non-lesional" = "darkolivegreen4", "Control" = "dodgerblue4")) +
        theme_classic() +
        theme(axis.title.x = element_blank()) +
        ylab(paste0(as.name(input$gene), " (CPM)")) +
        stat_compare_means(comparisons = list(c("Lesional", "Perilesional"), c("Perilesional", "Non-lesional"), c("Non-lesional", "Control"), c("Lesional", "Non-lesional"), c("Lesional", "Control"), c("Perilesional", "Control")))

      output$dot <- renderPlot({
        dot_plot
      })
    })

    if (debug) {
      message(paste("Dot Plot Time: ", dot_plot_time, "\n"))
    }

    # Dot plot download
    output$download_dot <- downloadHandler(
      filename = "DotPlot.pdf",
      content = function(file) {
        ggsave(filename = file, plot = dot_plot, width = 5 * length(input$gene), height = 5)
      }
    )

    # Correlation plot

    # doesn't work. find updated version
    # convert vsd to numeric matrix
    # cor_plot_time <- system.time({
    #   Cor_Plot <- heatmap(gene = input$gene, as.matrix(vsd)
    #
    #   output$Correlation <- renderPlot({
    #     Cor_Plot
    #   })
    # })

    # if (debug) {
    #   message(paste("Correlation Plot Time: ", cor_plot_time, "\n"))
    # }
  }

  # on submit
  observeEvent(input$submit, {
    req(data())
    message("submit triggered")
    if (valid_input(input$gene)) {
      generate_plots(input$gene)
    }
  })
}

# Run the application
shinyApp(ui = ui, server = server)
