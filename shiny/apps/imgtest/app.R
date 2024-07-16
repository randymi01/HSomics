library(shiny)
library(png)
library(shinyjs)

ui <- navbarPage(
    title = "Your App Title",
    useShinyjs(),
    mainPanel(
        img(src = "umap_plot.png", height = "50%", width = "100%", align = "right")
    )
)

server <- function(input, output) {
    # Add a class to the image
}

# Run the application
shinyApp(ui = ui, server = server)
