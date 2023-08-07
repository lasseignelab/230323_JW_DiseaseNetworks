library(shiny)
library(pheatmap)
library(grid)
library(stringr)
library(RColorBrewer)
library(circlize)
library(ggplot2)

# Load Data
load("data.Rdata") #object called data

# tf names
all_tfs <- rownames(data)

# Server
server <- function(input, output) {
  # Reactive values to store the selected tfs
  selected_tfs_data <- reactiveValues(selected_tfs = character(0))
  
  # Observe changes in the selected tfs
  observeEvent(input$tf_input_select, {
    selected_tfs_data$selected_tfs <- c("X", input$tf_input_select, "X")
  })
  
  # Function to create the heatmap plot
  generate_heatmap <- function(selected_tfs) {
    # Convert the comma-separated input to a vector of tf names
    selected_tfs <- unlist(strsplit(selected_tfs, "\\s*,\\s*"))
    
    # Subset the data based on selected tfs
    selected_data <- subset(data, rownames(data) %in% selected_tfs)
    
    # cleaning up names for plotting
    colnames(selected_data) <- gsub("_", " ", colnames(selected_data))
    colnames(selected_data) <- str_to_title(colnames(selected_data))
    
    # convert dataframe to matrix
    mat <- as.matrix(selected_data)
    
    # Plot heatmap
    pheatmap(
      mat,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      show_row_dend = FALSE,
      show_column_dend = FALSE,
      row_names_gp = gpar(fontface = "bold", fontfamily = "Helvetica"),
      row_title_gp = gpar(fontface = "bold"),
      show_row_names = TRUE,
      show_column_names = TRUE,
      column_names_gp = gpar(fontface = "bold", fontfamily = "Helvetica"),
      use_raster = TRUE,
      top_annotation = NULL
    )
  }
  
  # Event for the plot button
  observeEvent(input$plot_button, {
    selected_tfs <- paste(selected_tfs_data$selected_tfs, collapse = ", ")
    if (!is.null(selected_tfs) && nchar(selected_tfs) > 0) {
      output$heatmap_plot <- renderPlot({
        generate_heatmap(selected_tfs)
      })
    }
  })
}
