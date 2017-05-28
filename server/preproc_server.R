preproc_server <- function(input, output, session) {

  ns <- NS("preproc")
  ##########################
  # Reactive Values        #
  ##########################
  DB   <- reactiveValues(names=DB.ls(db))
  STUDY <- reactiveValues(action="", update=0, ori=NULL, preview=NULL,  
    expr=NULL, clinical=NULL, 
    studyName=NULL, species=NULL, genes=NULL,
    DE_p=NULL, DE_lfc=NULL, MCMC=NULL)
  SUMMARY <- reactiveValues(brief=data.frame(NULL), table=data.frame(NULL))

  ##########################
  # Validation             #
  ##########################
  validate.study <- function(study) {
    if(is.null(study$MCMC))  {
      stop(MSG.datasetInput.noinput)
    }
    studyName <- study$studyName
    if(is.null(studyName) || studyName == "") {
      stop(MSG.study.noname)
    }
    if(studyName %in% DB$names) stop(MSG.study.duplicate(studyName))  
  }

  ##########################
  # Observers              #
  ##########################
  # watch for tab change, get the newest list of all data
  observeEvent(input$tabChange, {DB$names <- DB.ls(db)}, 
               label="tab change")

  # watch for pvalue files upload
  observeEvent(input$pvalfile, {
    if (!is.null(input$pvalfile)) {
      STUDY$action <- STUDY.pval.upload
      STUDY$update <- STUDY$update + 1
      STUDY$species <- input$species
    }
  }, label="p-value file upload")

  # watch for expression file upload
  observeEvent(input$exprfile, {
    if (!is.null(input$exprfile)) {
      STUDY$action <- STUDY.expression.upload
      STUDY$update <- STUDY$update + 1
      STUDY$species <- input$species
    }
  }, label="expression file upload")

  # watch for clinical file upload
  observeEvent(input$clinical, {
    if (!is.null(input$clinical)) {
      STUDY$action <- STUDY.clinical.upload
      STUDY$update <- STUDY$update + 1
    }
  }, label="clinical file upload")

  # Setting study from file type
  observeEvent(STUDY$update, {
    if (STUDY$action == STUDY.pval.upload) {
      wait(session, "Parsing pvalue File...")
      inFile <- input$pvalfile
      tryCatch({
        pvalfile <- read.pData(inFile$datapath)
        STUDY$genes <- rownames(pvalfile)
        STUDY$DE_p <- pvalfile[,"pvalue"]
        STUDY$DE_lfc <- pvalfile[,"logFC"]
        brief <- data.frame(species=STUDY$species, NumGenes=length(STUDY$DE_p))
        SUMMARY$brief = brief
        SUMMARY$table <- cbind(STUDY$DE_p, STUDY$DE_lfc)
        rownames(SUMMARY$table) <- STUDY$genes
        colnames(SUMMARY$table) <- c("pvalue", "logFC")
      }, error=function(error) {
        sendErrorMessage(session, MSG.file.corrupted)
      })
      done(session)
    } else if (STUDY$action == STUDY.expression.upload) {
      wait(session, "Parsing Expression File...")
      inFile <- input$exprfile
      tryCatch({
        STUDY$expr <- read.rawData(inFile$datapath)
        STUDY$genes <- rownames(STUDY$expr)
      }, error=function(error) {
        sendErrorMessage(session, MSG.file.corrupted)
      })
      done(session)
    }  else if (STUDY$action == STUDY.clinical.upload) {
      wait(session, "Parsing Clinical File...")
      inFile <- input$clinical
      tryCatch({
        STUDY$clinical <- read.groupData(inFile$datapath)
        brief <- data.frame(species=STUDY$species, 
          NumSamples=ncol(STUDY$expr), NumGenes=nrow(STUDY$expr), 
          NumCase=sum(STUDY$clinical==1),
          NumControl=sum(STUDY$clinical==2))
        SUMMARY$brief = brief
      }, error=function(error) {
        sendErrorMessage(session, MSG.file.corrupted)
      })
      done(session)
    }
  }, label="setting STUDY$ori from file upload or selection")
  
  # DE analysis
  observeEvent(input$SingleDE, {
    wait(session, "DE analysis in single study, may take a while")
   try({
      SingleDERes <- indDE(data=STUDY$expr, 
        group=STUDY$clinical, data.type=input$dtype, case.label="2", 
        ctrl.label="1")
      DE_lfc <- SingleDERes[,"logFC"]
      DE_p <- SingleDERes[,"pvalue"]
      STUDY$DE_p <- DE_p
      STUDY$DE_lfc <- DE_lfc
      SUMMARY$table <- cbind(STUDY$DE_p, STUDY$DE_lfc)
      rownames(SUMMARY$table) <- STUDY$genes
      colnames(SUMMARY$table) <- c("pvalue", "logFC")
    }, session)
    done(session)
  })

  # Select DE genes
  observeEvent(input$SelectDE, {
    wait(session, "Select DE genes")
   try({
      DE_index <- deSelect(SUMMARY$table, p.cut=input$PCut, 
        topDE.number = input$DENum)
      STUDY$DE_index = DE_index
      print(db@dir)
    }, session)
    done(session)
  })

  # preprocess single studies
  observeEvent(input$preprocSingleStudy, {
    wait(session, "Preprocess single study, may take a while")
    try({
      set.seed(20170513) 
      #STUDY$MCMC <- bayesP(SUMMARY$table)
      # load results to save time
      load("/Users/lizhu/Box Sync/Genomics_new (LIZ86@pitt.edu)/CrossSpecies/DocFromCharles/package-0509/all_MCMCout.RData")
      if(input$species == "Human"){
        STUDY$MCMC <- hb_MCMCout
        }else{
          STUDY$MCMC <- mb_MCMCout
        }
    }, session)
    done(session)
  })

  # Save and Metadata
  observeEvent(input$saveStudy, {
    wait(session, "saving study")
    try({
      STUDY$studyName <- input$studyName
      #validate.study(STUDY)
      study_use <- new("Study1",
        studyName=STUDY$studyName,
        species=STUDY$species,
        DE_index=STUDY$DE_index,
        MCMC=STUDY$MCMC,
        matched_genes=rep(NA,nrow(STUDY$expr))
      )
      DB.save(db, study_use)
      sendSuccessMessage(session, paste("Study", study_use@studyName, "saved."))
      DB$names <- DB.ls(db)
      reset("species")
      reset("inputType")
      reset("pvalfile")
      reset("exprfile")
      reset("clinical")
      reset("studyName")
      reset("PCut")
      reset("DENum")
      SUMMARY$brief <- data.frame(NULL)
      SUMMARY$table <- data.frame(NULL)
    }, session)
    done(session)
  }, label="save study")

  ##########################
  # Render output/UI       #
  ##########################
  # brief
  output$brief <- DT::renderDataTable({
    table <- SUMMARY$brief
    if (is.null(table))
      "No file uploaded"
    else
      SUMMARY$brief
  })

  # summary
  output$summary <- DT::renderDataTable({
    table <- SUMMARY$table
    if (is.null(table))
      "No file uploaded"
    else
      SUMMARY$table
  })
 
}
