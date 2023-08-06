library(shiny)
library(pheatmap)
library(stringr)
library(RColorBrewer)
library(circlize)
library(ggplot2)

# Load Data
load("data.Rdata") # object called data

# tf names
all_tfs <- rownames(data)

# UI
ui <- fluidPage(
  # Add lab logo
  tags$img(src = "header.png", 
           style = "width: 100%; height: auto;"),
  
  # Information on Data and App -- Main Tab
  div(
    h3("Activity Score Calculation and Information"),
    p(
      "This application created by the",
      a(href = "https://www.lasseigne.org", "Lasseigne Lab", target = "_blank"), # opens in a new tab 
      "enables researchers to investigate the TF Activity of 758
      transcription factors (TFs) from the Genotype-Tissue Expression (GTEx) project
      across 31 human tissues.",
      br(),
      tags$figure(
        tags$img(
          src = "overview.png",
          width = "700px",
          height = "450px",
        )
      ),
      tags$ul(
        tags$li(
          "All scores were calculated using GTEx transcripts per million (TPM) expression data for 31 tissues (x-axis)."
        ),
        tags$li(
          "TPM data was combined with prior knowledge from",
          a(href = "https://github.com/saezlab/CollecTRI", "CollecTRI.", target = "_blank"),
        ),
        tags$li(
          "We used a multivariate linear model (run_mlm) with a minimum threshold of 5 targets per TF to calculate activity scores (represented as t-values) for all TFs",
          a(href = "https://saezlab.github.io/decoupleR/", "decoupleR.", target = "_blank")
        ),
        tags$li(
          "TF Activity scores will be represented on a scale from positive to negative corresponding to activation or inactivation of TFs.",
        ),
        tags$li(
          "We scaled and centered data prior to summarizing the averaged regulator activity for each TF.",         
        ),
        tags$li(
          "Positive and negative scores denote increased and decreased TF activity, respectively.",
        )
      ),
      p(
        "For more information about GTEx, visit their",
        a(href = "https://www.gtexportal.org/home/", "website", target = "_blank")
      )
    )
  ),  # Closing parenthesis added here to properly close the div

  div(
    p(
      strong("Instructions:"),
      br(),
      "1. Select one or more TFs from the list.",
      br(),
      "2. Click the 'Plot Heatmap' button to generate the visualization.",
      br(),
      em("Note: The full dataset of activity scores can be found",
         a(href = "https://github.com/lasseignelab/230323_JW_DiseaseNetworks/blob/main/results/tf_activity/tf_activity_scores.csv", "here", target = "_blank"))
    )
  ),

  # First fluid row for input elements
  fluidRow(
    column(
      width = 4, # Adjust the width to your preference for input elements
      # Input elements here (e.g., selectInput, actionButton, etc.)
      selectInput("tf_input_select", "Select TF(s):", choices = all_tfs, multiple = TRUE),
      actionButton("plot_button", "Plot Heatmap")
    )
  ),

  # Second fluid row for the plot
  fluidRow(
    column(
      width = 12, # The plot will occupy the full width of the page
      # Wrap the plotOutput in a div with a fixed height and scrollable content
      div(style = "overflow: auto;",
        # This column will contain the plot output
        plotOutput("heatmap_plot")
      )
    )
  )
)

