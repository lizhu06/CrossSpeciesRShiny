setting_ui <- function(id, label = "global settings") {
  ns <- NS(id)
  tabPanel("setting", value=id,
    mainPanel(
      h3("Session Information"),
      verbatimTextOutput(ns("urlText")),
      h3("Directory for Saving Output Files:", style="display:inline"),
      helpIcon("working_dir_help", HELP.working.dir),
      directoryInput(ns('directory'), label='select a directory')
    )
  )
}
