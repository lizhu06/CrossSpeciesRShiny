preproc_ui <- function(id, label= "preprocessing data") {
  ns <- NS(id)

  tabPanel("Upload", value=id,
    sidebarLayout(
      sidebarPanel(
        useShinyjs(),
        ##########################
        # Upload Data 
        ##########################
        #### chose species
        radioButtons(ns("species"), label = "Chose species",
                choices = list("Human" = "Human", "mouse" = "mouse"),
                selected = "Human"),

        #### chose input data type
        radioButtons(ns("inputType"), label = "Input data type",
                choices = list("p-value" = 1, "Raw data" = 2),
                selected = 1),

        #### if input data are p-values
        conditionalPanel(
          condition = paste0("input['", ns("inputType"), "'] == '1' "),
          fileInput(ns("pvalfile"), 'Upload p-values and logfc file',
            accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv'))
          ),

        #### if input data are raw data perform single DE first
        conditionalPanel(
          condition = paste0("input['", ns("inputType"), "'] == '2' "),
          radioButtons(ns("dtype"), label = "Data type",
                        choices = list("microarray" = "microarray", 
                          "RNAseq" = "RNAseq"),
                          selected = "microarray"),
          fileInput(ns("exprfile"), 'Upload raw data CSV file',
            accept=c('text/csv', 'text/comma-separated-values,text/plain', 
              '.csv')
          ),

          fileInput(ns("clinical"), 'Upload clinical CSV file',
            accept=c('text/csv', 'text/comma-separated-values,text/plain', 
              '.csv')
          ),
          
          actionButton(ns('SingleDE'), 'DE analysis', class="btn-success"),
          tags$hr(),

          p("Select DE genes:"),
          numericInput(ns("PCut"), label = "pvalue cut-off", 
           value = NA),
          numericInput(ns("DENum"), label = "Number of DE genes", 
            value = 1000),
          actionButton(ns('SelectDE'), 'Select DE genes', class="btn-success")
        ),
    
        #### calculate signed delta and q-values
        actionButton(ns('preprocSingleStudy'), 'preprocess single study', class="btn-success"),
        tags$hr(),

        ##########################
        # Save and Metadata      #
        ##########################
        #h4("Saving Study:"),
        textInput(ns("studyName"), "Save study:", "study name"),
        actionButton(ns('saveStudy'), 'save single study', icon=icon("save"), class="btn-success")
      ),


      mainPanel(
        h3("Study summary"),
        DT::dataTableOutput(ns("brief")),
        h3("DE summary"),
        DT::dataTableOutput(ns("summary"))
      )
    )
  )
}
