saved_data_ui <- function(id, label = "saved data of single study or multiple study") {
  ns <- NS(id)
  tabPanel("Saved Data", value=id,
    sidebarLayout(
      sidebarPanel(
        h4("Selected Datasets"), helpIcon(ns('merge_select'), HELP.select.datasets),
	      textOutput(ns("selected")),
        actionButton(ns('delete'), 'Delete Selected Data', icon=icon("trash"), 
          class="btn-danger"),
	      tags$hr(),
        fileInput(ns("orthologous"), 'Upload orthologous file',
          accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
        actionButton(ns("merge"), 'Match and merge', icon=icon("rocket"))
      ),
      mainPanel(
        h3("List of saved data"),
        DT::dataTableOutput(ns("table"))
      )
    )
  )
}
