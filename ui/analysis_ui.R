
analysis_ui <- function(id, label = "Analysis") {
  ns <- NS(id)
  tabPanel("Analysis", value=id,
    sidebarLayout(
      sidebarPanel(
        selectInput(ns("measure"), label = "Select measure", 
          choices = list("Youden" = "youden", "Fmeasure" = "Fmeasure", 
            "geo.mean" = "geo.mean"), 
          selected = "Fmeasure"),
        numericInput(ns("permNum"), 
          label = "Number of permutations for global ARS", 
          value = 100),
        actionButton(ns('globalARS'), 'Calculate global ARS', 
          class="btn-success"),
        tags$hr(),

        actionButton(ns('pathwayEnrich'), 
          'pathway enrichment analysis for each studies', 
          class="btn-success"),
        tags$hr(),

        numericInput(ns("pathwaysizeLowerCut"), 
          label = "Pathway size lower cut", value = 5),
        numericInput(ns("pathwaysizeUpperCut"), 
          label = "Pathway size upper cut", value = 199),
        numericInput(ns("overlapsizeCut"), 
          label = "Overlap size upper cut", value = 5),
        numericInput(ns("medDECut"), 
          label = "Median DE cut", value = 3),
        numericInput(ns("qfisherCut"), 
          label = "Fisher q-val cut", value = 0.05),
        numericInput(ns("B"), 
          label = "Number of permutations for pathway ARS", value = 100),
        actionButton(ns('pathwayARS'), 'Calculate pathway ARS', 
          class="btn-success"),
        tags$hr(),

        numericInput(ns("C"), 
          label = "Number of clusters", value = 7),
        actionButton(ns('pathwayClustering'), 'Pathway clustering', 
          class="btn-success"),
        tags$hr(),

        uiOutput(ns('selectPathwayName')),
        actionButton(ns('modelAnalysis'), 'Analyse the models', 
          class="btn-success"),
        tags$hr()

        ),
      mainPanel(
        h3("Global ARS is"),
        DT::dataTableOutput(ns("globalARSTable")),
        h3("Pathway enrichment results is"),
        DT::dataTableOutput(ns("pathwayEnrichTable")),
        h3("Pathway ARS results is"),
        DT::dataTableOutput(ns("pathwayARSTable")),
        h3("Pathway clustering MDS plot"),
        imageOutput(ns("MDSPlot"), height = 500),
        h3("Plot MDS plot for selected pathway"),
        imageOutput(ns("modelMDSPlot"), height = 500),
        h3("Model clustering for selected pathway"),
        imageOutput(ns("modelHeatmap"), height = 500),
        h3("Gene heatmap for selected pathway"),
        imageOutput(ns("geneHeatmap"), height = 500),
        h3("KEGG pathway viewer"),
        imageOutput(ns("KEGGtopo"), height = 500)
      )
    )
  )
}
