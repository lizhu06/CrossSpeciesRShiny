analysis_server <- function(input, output, session){
  ns <- NS("analysis")
  
  ##########################
  # Reactive Values        #
  ##########################
  DB <- reactiveValues(MergedDB=MergedDB.load(db), 
    dat=NULL, deIndex=NULL, 
    globalARI=NULL, globalARSpvalue=NULL, 
    metaSummary=data.frame(NULL), pathARS_summary=data.frame(NULL),
    pathARS=NULL, pathARSpvalue=NULL, pathway.list=NULL
    )

  ##########################
  # Observers              #
  ##########################
  # tab change to load merged data
  observeEvent(input$tabChange, {
    DB$MergedDB <- MergedDB.load(db)
    print(DB.load.working.dir(db))
    print(length(DB$MergedDB))
  })

  # calculate global ARS
  observeEvent(input$globalARS, {
    if (length(DB$MergedDB) > 1) {
      wait(session, "Calculate global ARS, may take a while")

      # organize the data list
      all_species <- sapply(1:length(DB$MergedDB), function(x) 
        DB$MergedDB[[x]]@species)

      DB$dat <- de <- DB$deIndex <- list()
      de <- list()
      for(i in 1:length(DB$MergedDB)){
        genes <- DB$MergedDB[[i]]@matched_genes
        DB$dat[[i]] <- (DB$MergedDB[[i]]@MCMC)[genes,]
        de[[i]] <- intersect(names(DB$MergedDB[[i]]@DE_index), 
          rownames(DB$dat[[i]]))
        DB$deIndex[[i]] <- which(rownames(DB$dat[[i]])%in%de[[i]])
      }

      # global ARS (directly show results to save time)
      #DB$globalARS <- ars(DB$dat[[1]], DB$dat[[2]], DB$deIndex[[1]], 
      #  DB$deIndex[[2]], measure=input$measure)
      #globalARSperm <- permGlobal(DB$dat[[1]], DB$dat[[2]], B=input$permNum)
      #DB$globalARSpvalue <- arsPvalue(DB$globalARS, globalARSperm)
      load("/Users/lizhu/Box Sync/Genomics_new (LIZ86@pitt.edu)/CrossSpecies/DocFromCharles/package-0509/hbmb_globalARS.RData")
      DB$globalARS <- hbmb_globalARS
      load("/Users/lizhu/Box Sync/Genomics_new (LIZ86@pitt.edu)/CrossSpecies/DocFromCharles/package-0509/hbmb_globalARSpvalue.RData")
      DB$globalARSpvalue <- hbmb_globalARSpvalue
      
      done(session)
    }
  })

  # pathway enrichment analysis
  observeEvent(input$pathwayEnrich, {
    if(length(DB$MergedDB) > 1) {
      wait(session, "Pathway enrichment")
      
      # pathway enrichment analysis
      DB$pathway.list <- c(GOBP.genesets, GOCC.genesets, GOMF.genesets,
                        Reactome.genesets, Immunologic.genesets,
                        KEGG.genesets, Biocarta.genesets) 
      rownames(DB$dat[[2]]) <- rownames(DB$dat[[1]])

      #(directly load results to save time!!!!!!!!)
      #pathSummary <- list()
      #for(i in 1:length(DB$MergedDB)){
      #  pathSummary[[i]] <- pathEnrich(dat=DB$dat[[i]], 
      #    deIndex=DB$deIndex[[i]],
      #    pathway.list=pathway.list)
      #}
      #DB$metaSummary <- metaPath(pathSummary=pathSummary, 
      #  pathway.list=pathway.list)

      load("/Users/lizhu/Box Sync/Genomics_new (LIZ86@pitt.edu)/CrossSpecies/DocFromCharles/package-0509/pathSummary.RData")
      load("/Users/lizhu/Box Sync/Genomics_new (LIZ86@pitt.edu)/CrossSpecies/DocFromCharles/package-0509/metaSummary.RData")
      DB$metaSummary <- metaSummary

      done(session)
    }
  })

  # pathway ARS
  observeEvent(input$pathwayARS, {
    print(length(DB$MergedDB))
    if(length(DB$MergedDB) > 1) {
      wait(session, "Pathway enrichment analysis")
      
      selectPathways <- pathSelect(DB$metaSummary, 
        input$pathwaysizeLowerCut, input$pathwaysizeUpperCut, 
        input$overlapsizeCut, input$medDECut, input$qfisherCut)
      K <- length(selectPathways)

      DB$pathARS <- rep(NA,K)
      names(DB$pathARS) <- selectPathways
      dat.gene <- rownames(DB$dat[[1]]) ## same for all data (already matched)
      for(k in 1:K){
        print(k)
        path.gene <- intersect(DB$pathway.list[[selectPathways[k]]], dat.gene)
        path.dat <- lapply(1:length(DB$dat), function(x) DB$dat[[x]][path.gene, ])
        de.gene <- lapply(1:length(DB$dat), function(x) dat.gene[DB$deIndex[[x]]])
        
        path.deIndex <- list()
        for(i in 1:length(DB$dat)){
          if(length(intersect(de.gene[[i]] , path.gene))<=3 ){
            path.deIndex[[i]] <- 1:nrow(path.dat[[i]])
          }else{
            path.deIndex[[i]] <- which(rownames(path.dat[[i]])%in% 
              intersect(de.gene[[i]], path.gene))
          }   
        }

        DB$pathARS[k] <- ars(path.dat[[1]], path.dat[[2]], 
          path.deIndex[[1]], path.deIndex[[2]], measure="Fmeasure")
      }

      ## pathway ARS permutation
      #pathARSperm <- matrix(NA, nrow=input$B, ncol=K)
      #colnames(pathARSperm) <- selectPathways
      #pathwaySize <- sapply(1:K, FUN=function(k){
      #  length(intersect(rownames(dat[[1]]), 
      #    pathway.list[[selectPathways[k]]]))
      #})
      #pathARSperm <- permPathway(DB$dat[[1]], DB$dat[[2]], selectPathways, pathwaySize, B=input$B)

      #DB$pathARSpvalue <- sapply(1:K,function(k) 
      #  {arsPvalue(DB$pathARS[k], pathARSperm[,k])})
      #names(DB$pathARSpvalue) <- selectPathways
      #DB$pathARS_summary <- cbind(DB$pathARS, DB$pathARSpvalue)

      # directly load results to save time
      load("/Users/lizhu/Box Sync/Genomics_new (LIZ86@pitt.edu)/CrossSpecies/DocFromCharles/package-0509/hbmb_pathARS_summary.RData")
      DB$pathARS_summary <- hbmb_pathARS_summary

      # pathway ARS table
      output$pathwayARSTable <- DT::renderDataTable({
        if (!(is.null(DB$pathARS_summary))) {
          DB$pathARS_summary
        }   
      })

      done(session)
    }
  })


  # pathway clustering
  observeEvent(input$pathwayClustering, {
    if(length(DB$MergedDB) > 1) {
      wait(session, "Pathway clustering")
      
      results <- clustPathway(ARSFpvalue)
      cluster.assign <- results[[input$C]]$consensusClass
      scatter.index <- scatter(-log10(ARSFpvalue), cluster.assign)
      mds.res <- mdsPathway(arsPvalue=ARSFpvalue,
                           cluster.assign=cluster.assign,
                           scatter.index=scatter.index, 
                           folder=DB.load.working.dir(db))
      output$MDSPlot <- renderImage({
        img.src <- paste(DB.load.working.dir(db), 
          "/MDS_pathway_", input$C, "_clusters.png", sep="")
        list(src=img.src, contentType='image/png', alt="module")
      }, deleteFile = FALSE)
      output$selectPathwayName = renderUI({
          selectInput(ns('selectPathwayName'), 'select pathway to analyze', 
            rownames(ARSFpathway_out),
           selected=rownames(ARSFpathway_out)[1])
      })
      done(session)
    }
  })

  # Per-pathway model analysis
  observeEvent(input$modelAnalysis, {
    if(length(DB$MergedDB) > 1) {
      wait(session, "Model analysis")
      
      ### MDS plot of models for selected pathways
      model.name <- c("ht","hb","hi","hs","hl","ha",
                      "mt","mb","mi","ms","ml","ma")
      arsPath <- ARSFpathway_out[input$selectPathwayName,] ## the name of pairs has to be in the form: paste(name1,name2,sep="")
      res <- mdsModel(arsPath,model.name,input$selectPathwayName)
      png(paste(DB.load.working.dir(db), "/mds_", input$selectPathwayName,".png", sep=""))
      print(res)
      dev.off()

      output$modelMDSPlot <- renderImage({
        img.src <- paste(DB.load.working.dir(db), 
          "/mds_", input$selectPathwayName, ".png", sep="")
        list(src=img.src, contentType='image/png', alt="module")
      }, deleteFile = FALSE)

      ### clustering heatmap
      arsPvalue <- ARSFpvalue[input$selectPathwayName,] 
      path.cluster.assign <- clustModel(arsPvalue,model.name)
      res <- heatmapModel(arsPvalue,model.name,path.cluster.assign,
                          input$selectPathwayName, 
                          folder=DB.load.working.dir(db))
      output$modelHeatmap <- renderImage({
        img.src <- paste(DB.load.working.dir(db), "/", input$selectPathwayName,'_clustHeatmap.png',sep="")
        list(src=img.src, contentType='image/png', alt="module")
      }, deleteFile = FALSE)

      ## gene heatmap
      signPM.list <- list(apply(DB$dat[[1]],1,mean),apply(DB$dat[[2]],1,mean))
      names(signPM.list) <- c("hb","mb")
      res <- heatmapGene(signPM.list, pathway.list=DB$pathway.list, 
                         pathway.name=input$selectPathwayName, 
                         folder=DB.load.working.dir(db))
      output$geneHeatmap <- renderImage({
        img.src <- paste(DB.load.working.dir(db), '/Heatmap', 
          input$selectPathwayName,'.png', sep="")
        list(src=img.src, contentType='image/png', alt="module")
      }, deleteFile = FALSE)

      ### KEGG viewer
      if(grepl("KEGG",input$selectPathwayName) == TRUE){
        overlap.genes <- intersect(rownames(DB$dat[[1]]), 
          DB$pathway.list[[input$selectPathwayName]])
        signPM.mat <- cbind(apply(DB$dat[[1]][overlap.genes,],1,mean),
                            apply(DB$dat[[2]][overlap.genes,],1,mean))
        colnames(signPM.mat) <- c("hb","mb")
        res <- keggView(mat=signPM.mat,kegg.name=input$selectPathwayName, folder=DB.load.working.dir(db))

        files <- list.files(paste(DB.load.working.dir(db), "/topo_", 
            input$selectPathwayName, sep=""))
        target_file <- files[grepl("multi.png", files)]
        output$KEGGtopo <- renderImage({
          img.src <- paste(DB.load.working.dir(db), "/topo_", 
            input$selectPathwayName, "/", target_file,sep="")
          list(src=img.src, contentType='image/png', alt="module")
        }, deleteFile = FALSE)

        print(paste(DB.load.working.dir(db), "/topo_", 
            input$selectPathwayName, "/", target_file,sep=""))
      }
      

      done(session)
    }
  })

  ##########################
  # Render output/UI       #
  ##########################
  # global ARS
  output$globalARSTable <- DT::renderDataTable({
    #table <- summary(DB$final_data[[1]]@MCMC)
    if (is.null(DB$globalARS)){
        data.frame(NULL) 
    }else{
      res <- matrix(1, length(DB$MergedDB), 
        length(DB$MergedDB))
      colnames(res) <- rownames(res) <-sapply(1:length(DB$MergedDB), 
        function(x) DB$MergedDB[[x]]@studyName)
      res[1,2] <- res[2,1] <- paste(round(DB$globalARS,2), " (p-val= ", 
        round(DB$globalARSpvalue,2), ")", sep="")
      res
    }
      
  })

  # pathway enrichment analysis
  output$pathwayEnrichTable <- DT::renderDataTable({
    if (!(is.null(DB$metaSummary))) {
      DB$metaSummary
    }   
  })



}