
saved_data_server <- function(input, output, session) {

  ns <- NS("saved_data")
 
  DB <- reactiveValues(meta=meta(db), 
    all_studies=DB.load(db, list.files(path=db@dir)))

  observeEvent(input$tabChange, {
    DB$meta <- meta(db)
    DB$all_studies <- DB.load(db, list.files(path=db@dir))
  })

  output$table <- DT::renderDataTable(DT::datatable({
    DB$meta
  }))

  output$selected <- renderText({
    selected <- input$table_rows_selected
    if(length(selected) == 0)
      "You haven't select anything yet"
    else
      paste(rownames(meta(db)[selected,]), sep=", ")
    
  })

  observeEvent(input$delete, {
    selected <- input$table_rows_selected
    if(length(selected) == 0) {
      sendErrorMessage(session, "You haven't select anything yet")
    } else {
      selected <- rownames(meta(db)[selected,])
      DB.delete(db, selected)
      sendSuccessMessage(session, paste(selected, "deleted"))
      DB$meta <- meta(db)
      DB$all_studies <- DB.load(db, list.files(path=db@dir))
    }
  })

  observeEvent(input$orthologous, {
    if (!is.null(input$orthologous)) {
      inFile <- input$orthologous
      DB$full_ortholog <- read.csv(inFile$datapath, stringsAsFactors = F,
         header=T)
    }
  })

  observeEvent(input$merge, {
    wait(session, "Match and merge")
    try({
        print(length(DB$all_studies)) 
        if(nrow(DB$meta)==2 && length(unique(DB$meta[,"species"]))==2){
          out <- upperlowerMatch(rownames(DB$all_studies[[1]]@MCMC), 
            rownames(DB$all_studies[[2]]@MCMC))
          out2 <- orthMatch(out[[4]],out[[3]],
                            ortholog=DB$full_ortholog, reference="2")
          DB$all_studies[[1]]@matched_genes <- c(out[[1]], out2[[2]])
          DB$all_studies[[2]]@matched_genes <- c(out[[2]], out2[[1]])
          MergedDB <- DB$all_studies
          saveRDS(MergedDB, file=paste(DB.load.working.dir(db), "MergedDB.rds", sep="/"))
        }
        message = paste("Data are successfully merged")
        sendSuccessMessage(session, message)
    })
    done(session)
  }, label="save study")

}
