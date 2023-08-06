library(shiny)
library(pheatmap)
library(stringr)
library(RColorBrewer)
library(circlize)
library(ggplot2)

source("ui.R")
source("server.R")

# Run the app
shinyApp(ui, server)