pkg_time <- system.time({
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
  library(shinyjs)
})

ui <- navbarPage(
  useShinyjs(),
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
    p("GSE154775, GSE155850"),
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
          h3("Select Plots"),
          checkboxInput(
            inputId = "feature_plot",
            label = strong("Feature Plot"),
            value = TRUE
          ),
          checkboxInput(
            inputId = "vln_plot",
            label = strong("Violin Plot"),
            value = TRUE
          ),
          checkboxInput(
            inputId = "dot_plot",
            label = strong("Dot Plot"),
            value = TRUE
          ),
          checkboxInput(
            inputId = "cor_plot",
            label = HTML("<strong>Correlation Plot</strong> <span style='color:red;'> (~1 min) </span>"),
            value = FALSE
          ),
          actionButton("submit", "Submit")
        ),
        mainPanel(
          fluidRow(
            div(
              style = "overflow: hidden;",
              imageOutput("dr_plot", width = "100%")
            )
          ),
          fluidRow(
            column(6, plotOutput("feature")),
            column(6, plotOutput("vln"))
          ),
          fluidRow(
            column(6, downloadLink("download_feature", "Download PDF")),
            column(6, downloadLink("download_vln", "Download PDF"))
          ),
          br(),
          fluidRow(
            column(12, plotOutput("dot"))
          ),
          fluidRow(
            column(12, downloadLink("download_dot", "Download PDF"))
          ),
          fluidRow(
            column(12, plotOutput("Correlation"))
          ),
          fluidRow(
            column(12, downloadLink("download_cor", "Download PDF"))
          )
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
  tags$footer(
    tags$div(
      style = "background-color: #f5f5f5; padding: 5px; text-align: center; font-size: 10px; border-top: 1px solid #ddd; position: fixed; left: 0; bottom: 0; width: 100%;",
      "This portal is developed and maintained by the Center for Bioinformatics, Department of Public Health Sciences and the Center for Cutaneous Biology and Immunology Research at Henry Ford Health, Detroit, Michigan."
    )
  )
)

heatmap <- function(gene = "NULL", vsd) {
  gene <- as.character(gene)
  counts <- t(assay(vsd))

  # first selects all columns except gene of interest, then runs a correlation test
  # for each non-interest genes and the gene of interest
  b <- apply(counts[, !(colnames(counts) == gene)], 2, function(x) {
    cor.test(counts[, gene], x)
  })
  cor.vals <- sapply(b, "[[", "estimate")
  p.vals <- sapply(b, "[[", "p.value")
  FDR <- p.adjust(p.vals, method = "bonferroni")
  tmp <- as.data.frame(cbind(p.vals, cor.vals, FDR))
  tmp <- tmp[which(tmp$FDR <= 0.05 & abs(tmp$cor.vals) >= 0.5), ]
  tissue_type <- as.data.frame(colData(vsd)[, c("Run", "tissue_type")])
  counts1 <- t(counts)[which(colnames(counts) %in% rownames(tmp)), order(assay(vsd)[colnames(counts) %in% gene, ], decreasing = T)]
  tissue_type <- tissue_type[match(colnames(counts1), tissue_type$Run), c(2)]
  tissue_type <- as.data.frame(tissue_type)
  rownames(tissue_type) <- colnames(counts1)
  pheatmap(counts1,
    cluster_rows = TRUE, show_rownames = FALSE, show_colnames = F,
    cluster_cols = FALSE, annotation_col = tissue_type, main = gene
  )
}

server <- function(input, output, clientData) {
  # reactive will dynamically reload data when RDS is updated
  # no need since will be refreshed each time.

  # is reactive required for the data?

  withProgress(message = "Please Wait, Loading Data", value = 0, {
    data <- list()

    incProgress(0.1, detail = "Reading Seurat Objects...")

    data$seurat <- readRDS("../../data/hsomicsdb/seurat/HS_merged_lite.RDS")
    incProgress(0.3, detail = "Reading Adjusted CPM...")

    # bulk
    data$cp <- readRDS("../../data/hsomicsdb/Combat_Adjusted_CPM.rds")
    incProgress(0.3, detail = "Reading Markers...")

    data$mk <- read.delim("../../data/hsomicsdb/markers/hs.markers.top30_genes.csv",
      sep = ",", stringsAsFactors = FALSE
    )
    incProgress(0.3, detail = "Reading VSD...")

    # bulk
    data$vsd <- readRDS("../../data/hsomicsdb/Updated_VSD.rds")

    incProgress(0.1, detail = "Done")
  })


  # check if gene exists in seurat object
  valid_input <- function(gene) {
    if (gene %in% rownames(data$seurat)) {
      return(TRUE)
      output$warning_message <- renderText({
        " "
      })
    } else {
      output$warning_message <- renderText({
        paste("Gene ", gene, " not found in the dataset")
      })
      return(FALSE)
    }
  }


  # return pregenerated UMAP
  # reactive container for when width changes
  observe({
    req(data)
    output$dr_plot <- renderImage(
      {
        list(
          src = "umap_plot.png",
          width = "100%",
          height = "100%"
        )
      },
      deleteFile = FALSE
    )
  })

  # hide download buttons on startup:
  shinyjs::hide("download_feature")
  shinyjs::hide("download_vln")
  shinyjs::hide("download_dot")
  shinyjs::hide("download_cor")

  clear_plots <- function() {
    output$feature <- renderPlot({})
    output$vln <- renderPlot({})
    output$dot <- renderPlot({})
    output$Correlation <- renderPlot({})
    if (input$feature_plot) {
      shinyjs::show("download_feature")
    } else {
      shinyjs::hide("download_feature")
    }

    if (input$vln_plot) {
      shinyjs::show("download_vln")
    } else {
      shinyjs::hide("download_vln")
    }

    if (input$dot_plot) {
      shinyjs::show("download_dot")
    } else {
      shinyjs::hide("download_dot")
    }

    if (input$cor_plot) {
      shinyjs::show("download_cor")
    } else {
      shinyjs::hide("download_cor")
    }
  }

  generate_plots <- function(gene) {
    # make sure datasets have been loaded
    req(data)

    # Feature plot
    if (input$feature_plot) {
      feature_plot <- Seurat::FeaturePlot(data$seurat,
        features = input$gene,
        reduction = "umap",
        pt.size = 1, raster = FALSE, cols = c(viridis(5)),
        order = T
      )
      output$feature <- renderPlot({
        req(input$submit)
        feature_plot
      })
      incProgress(0.2, detail = "Processing...")

      # Feature plot download
      output$download_feature <- downloadHandler(
        filename = "feature.pdf",
        content = function(file) {
          ggsave(filename = file, plot = feature_plot(), width = 5 * length(input$gene), height = 4)
        }
      )
    }

    # Violin plot
    if (input$vln_plot) {
      vln_plot <- Seurat::VlnPlot(data$seurat,
        features = input$gene,
        pt.size = input$point_size
      ) + NoLegend() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

      output$vln <- renderPlot({
        req(input$submit)
        vln_plot
      })
      incProgress(0.2, detail = "Processing...")
      # Violin plot download
      output$download_vln <- downloadHandler(
        filename = "violin.pdf",
        content = function(file) {
          ggsave(filename = file, plot = vln_plot(), width = 7 * length(input$gene), height = 5)
        }
      )
    }

    # Dot plot
    if (input$dot_plot) {
      dot_plot <- ggplot(data$cp, aes(
        x = factor(tissue_type, levels = c("Lesional", "Perilesional", "Non-lesional", "Control")),
        y = !!rlang::sym(input$gene),
        fill = tissue_type
      )) +
        geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5, stroke = NA) +
        scale_fill_manual(values = c(
          "Lesional" = "firebrick4",
          "Perilesional" = "darkgoldenrod3",
          "Non-lesional" = "darkolivegreen4",
          "Control" = "dodgerblue4"
        )) +
        theme_classic() +
        theme(axis.title.x = element_blank()) +
        ylab(paste0(input$gene, " (CPM)")) +
        stat_compare_means(comparisons = list(
          c("Lesional", "Perilesional"),
          c("Perilesional", "Non-lesional"),
          c("Non-lesional", "Control"),
          c("Lesional", "Non-lesional"),
          c("Lesional", "Control"),
          c("Perilesional", "Control")
        )) +
        ggtitle(input$gene) +
        theme(
          plot.title = element_text(
            size = 16,
            face = "bold",
            hjust = 0.5
          ),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 10)
        )

      output$dot <- renderPlot({
        req(input$submit)
        dot_plot
      })
      incProgress(0.2, detail = "Processing...")

      # Dot plot download
      output$download_dot <- downloadHandler(
        filename = "DotPlot.pdf",
        content = function(file) {
          ggsave(filename = file, plot = dot_plot, width = 5 * length(input$gene), height = 5)
        }
      )
    }

    # correlation plot
    if (input$cor_plot) {
      # takes 30ish seconds
      cor_plot <- heatmap(gene = input$gene, data$vsd)
      output$Correlation <- renderPlot({
        cor_plot
      })
      output$download_cor <- downloadHandler(
        filename = "CorPlot.pdf",
        content = function(file) {
          ggsave(filename = file, plot = cor_plot, width = 5 * length(input$gene), height = 5)
        }
      )
    }
  }



  # on submit
  observeEvent(input$submit, {
    req(data)
    if (valid_input(input$gene)) {
      clear_plots()
      withProgress(message = "Generating plots", value = 0, {
        incProgress(0.1, detail = "Processing...")
        generate_plots(input$gene)
        incProgress(0.5, detail = "Done")
      })
    }
  })
}

# Run the application
shinyApp(ui = ui, server = server)
