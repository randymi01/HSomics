library(Seurat)
library(shiny)
library(DT)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(pals)
library(dplyr)
library(bslib)
library(DESeq2)
library(pheatmap)

ui <- navbarPage(title = "HS-OmicsDB: Hidradenitis Suppurativa Omics Database",
                 
                 # Output plots
                 tabPanel(title = "Home",
                          br(),
                          hr(),
                          #h4(strong("Project Description")),
                          p(style="text-align: justify; font-size = 25px",
                            "Here we present HS 'omics database (HS-OmicsDB)', 
          an initiative to integrate the growing body of 'omics data' 
          collected from studies of HS. We hope this initiative provides
          a tool for the HS research community to enhance investigations 
          into the etiology of and specific functions of novel disease driving genes. 
          Our database compiles publicly available high-throughput bulk-RNA and single 
          cell RNA sequencing studies with publicly available fastq files from studies 
          on patients with HS deposited on the National institute of Health gene expression omnibus (GEO)."),
                          p("Users can explore genes of interest to determine tissue-specific and cell-specific gene expression: "), 
p("Current data from bulk-sequencing studies includes: "),
tags$ul(
tags$li("73 lesional (L)"),
tags$li("28 perilesional (PL)"),
tags$li("30 nonlesional (NL)"),
tags$li("19 healthy control (HC)")),                          
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
                 
                 tabPanel(title = "Visualize gene expression",
                          fluidPage(
                            sidebarLayout(
                              sidebarPanel(width = 3,
                                           h3("Data"),
                                           
                                           textInput("gene", "Genes", 'IL17A'),
                                           h3("Violin plot"),
                                           numericInput("pt_size", "Point size", value = 0.1, min = 0, max = 1, step = 0.1)),
                            
                            mainPanel(plotOutput("dr_plot", width = "100%"),
                                      fluidRow(
                                        splitLayout(cellWidths = c("50%", "50%"), plotOutput("feature"), plotOutput("vln"))
                                      ),
                                      fluidRow(
                                        splitLayout(cellWidths = c("50%", "50%"), downloadLink("download_feature", "Download PDF"), downloadLink("download_vln", "Download PDF"))
                                      ),
                                      plotOutput("dot"),
                                      downloadLink("download_dot", "Download PDF"),
                                      #plotOutput("vln"),
                                      #downloadLink("download_vln", "Download PDF"),
				      plotOutput("cor"),
				      downloadLink("download_cor","Download PDF")
                            )
                                       )),
                 ),        
                 tabPanel(title = "Cluster info",
                          tableOutput("ncell_table")),
                 
                 tabPanel(title = "Cluster markers",
                          ##p("Use the navigation tools to search, filter, and order the table of gene markers."),
                          DT::dataTableOutput("markers")),
                 
                 
                 h6("This portal is developed and maintained by the Center for Bioinformatics, Department of Public Health Sciences and the Center for Cutaneous Biology and Immunology Research at Henry Ford Health, Detroit, Michigan.")
)

server <- function(input, output, session) {
    heatmap <- function(gene = 'NULL', vsd){
  gene <- as.character(gene)
  counts <- t(assay(vsd))
  b = apply(counts[,!colnames(counts) %in% gene],2, function(x) {
    cor.test(counts[,gene],x)})
  cor.vals <- sapply(b, "[[", "estimate")
  p.vals <- sapply(b, "[[", "p.value")
  FDR = p.adjust(p.vals, method = 'bonferroni')
  FDR <- p.adjust(p.vals, method = 'bonferroni')
  tmp = as.data.frame(cbind(p.vals, cor.vals, FDR))
  tmp <- tmp[which(tmp$FDR <= 0.05 & abs(tmp$cor.vals) >= 0.5) ,]
  tissue_type <- as.data.frame(colData(vsd)[,c("Run","tissue_type")])
  counts1 <- t(counts)[which(colnames(counts) %in% rownames(tmp)),order(assay(vsd)[colnames(counts) %in% gene,], decreasing = T)]
  tissue_type <- tissue_type[match(colnames(counts1),tissue_type$Run),c(2)]
  tissue_type <- as.data.frame(tissue_type)
  rownames(tissue_type) <- colnames(counts1)
  pheatmap(counts1, cluster_rows=TRUE, show_rownames=FALSE, show_colnames = F,
           cluster_cols=FALSE, annotation_col=tissue_type) 
}
    #observe({
    #    
    #    sample <- input$sample
    #    
    #    if (is.null(sample)) sample <- character(0)
    #    
    #    updateSelectInput(session, "gene",
    #                     choices = genes[[sample]])
    #    
    #})
    
    # Reactive expressions
    #gene   <- reactive({input$gene})
    seurat <- reactive({readRDS("/shiny/data/hsomicsdb/seurat/HS_merged_lite.RDS")})
    cp <- reactive({readRDS("/shiny/data/hsomicsdb/Combat_Adjusted_CPM.rds")})
    mk     <- reactive({read.delim(paste0("/shiny/data/hsomicsdb/markers/hs.markers.top30_genes.csv"), sep = ',',
                                   stringsAsFactors = FALSE)})
    vsd <- reactive({readRDS('/shiny/data/hsomicsdb/Updated_VSD.rds')})
    
    # UMAP
    output$dr_plot <- renderPlot({Seurat::DimPlot(seurat(), reduction = 'umap',
                                                  pt.size = 1.1, label = T, cols = as.character(polychrome(34)),
                                                  raster = FALSE, label.size = 6) + NoLegend()})#, height = 800, width = 1200)
    
    # Number of cells
    output$ncell_table <- renderTable({
        n_cells <- as.data.frame(table(seurat()@active.ident))
        colnames(n_cells) <- c("Cluster", "Number of cells")
        return(n_cells)
    })
    
    # Markers
    output$markers <- DT::renderDataTable({
        mk() %>%
            dplyr::select(cluster, gene, avg_log2FC, p_val_adj, pct.1, pct.2) %>% 
            DT::datatable(filter = "top",
                          rownames = FALSE,
                          selection = "none") %>% 
            formatStyle("cluster",  fontWeight = "bold") %>% 
            formatStyle("gene",  fontWeight = "bold")
    })
    
    # Feature plot
    feature_plot <- reactive({Seurat::FeaturePlot(seurat(), features = input$gene,
                                               reduction = 'umap',
                                               pt.size = 1, raster = FALSE,cols = c(viridis(5)))})
    
    output$feature <- renderPlot({feature_plot()})#, height = 400, width = 600)
    
    # Feature plot download
    output$download_feature <- downloadHandler(filename = "feature.pdf",
                                               content = function(file) {
                                                   ggsave(filename = file, plot = feature_plot(), width = 5*length(input$gene), height = 4)
                                               })
    
    # Violin plot
    vln_plot <- reactive({Seurat::VlnPlot(seurat(), features = input$gene,
                                          pt.size = input$point_size) + NoLegend() +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))})
    
    output$vln <- renderPlot({vln_plot()})#, height = 400, width = 600)
    
    # Violin plot download
    output$download_vln <- downloadHandler(filename = "violin.pdf",
                                           content = function(file) {
                                               ggsave(filename = file, plot = vln_plot(), width = 7*length(input$gene), height = 5)
                                           })
    dot_plot <- reactive({ggplot(cp(), aes(x=factor(tissue_type, levels=c('Lesional','Perilesional','Non-lesional','Control')), y=eval(as.name(input$gene)), fill = tissue_type)) + 
                                           geom_dotplot(binaxis='y', stackdir='center',dotsize = 0.5,stroke=NA) + scale_fill_manual(values=c('Lesional' = "firebrick4", 'Perilesional' = "darkgoldenrod3", 'Non-lesional' = "darkolivegreen4",'Control' = "dodgerblue4")) +
                                           theme_classic()+ theme(axis.title.x = element_blank())+ ylab(paste0(as.name(input$gene), " (CPM)")) +
											stat_compare_means(comparisons = list( c('Lesional','Perilesional'), c('Perilesional','Non-lesional'), c('Non-lesional','Control'),c('Lesional','Non-lesional'), c('Lesional','Control'), c('Perilesional','Control') ))})    
    output$dot <- renderPlot({dot_plot()})
    
    output$download_dot <- downloadHandler(filename = "DotPlot.pdf",
                                           content = function(file) {
                                               ggsave(filename = file, plot = dot_plot, width = 5*length(input$gene), height = 5)
                                           })

	Cor_Plot <- reactive({heatmap(gene = input$gene, vsd)})
        output$cor <- renderPlot({Cor_Plot()})

    output$download_cor <- downloadHandler(filename = "CorPlot.pdf",
                                          content = function(file) {
                                             ggsave(filename = file, plot = dot_plot, width = 5*length(input$gene), height = 5)
                                        })

    
}

# Run the application 
shinyApp(ui = ui, server = server)