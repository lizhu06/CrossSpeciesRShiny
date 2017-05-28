# This file will be executed prior to app startup to setup the necessary environment
GLOBAL.network <- T
tryCatch({
  source("https://bioconductor.org/biocLite.R")
}, error=function(error){
  GLOBAL.network <<- F
})

# Try again with http if https fail
if (GLOBAL.network == F) {
  tryCatch({
    source("http://bioconductor.org/biocLite.R")
    GLOBAL.network <<- T
  }, error=function(error){
    GLOBAL.network <<- F
  })
}

if(GLOBAL.network) {
  installed <- installed.packages()[,"Package"]
  for (package in c("utils", "DMwR", "devtools", "DT", "shinyBS", 
    "Rcpp", "RcppArmadillo", "RcppGSL", "snowfall", "cvTools",
    "samr", "MASS", "ggrepel", "ggplot2", "gplots", "WGCNA",
    "tightClust")) {
    if (!(package %in% installed)) {
      install.packages(package, repos='http://cran.us.r-project.org')
    }
  }
  if (!("AnnotationDbi" %in% installed)) {
    biocLite("AnnotationDbi")
  }
  if (!("limma" %in% installed)) {
    biocLite("limma")
  }
  if (!("preproc" %in% installed)) {
    devtools::install_github("metaOmic/preproc")
  }
}

library(preproc)
library(shiny)
library(shinyBS)

library(Rcpp)
library(RcppArmadillo)
library(RcppGSL)
library(snowfall)
library(cvTools)
library(samr)
library(limma) # biocon
library(MASS)
library(ggrepel)
library(ggplot2)
library(gplots)
library(WGCNA)
library(tightClust)
library(ConsensusClusterPlus)
library(biomaRt)
library(KEGG.db)
library(org.Hs.eg.db)
library(pathview)

data(preproc.option) ##???

source("src/internal.R")
sourceCpp("src/BayesP.cpp")
sourceCpp("src/ARS.cpp")
source("src/readData.R")
source("src/indDE.R")
source("src/bayesP.R")
source("src/upperlowerMatch.R")
source("src/orthMatch.R")
source("src/ars.R")
source("src/perm.R")
source("src/arsPvalue.R")
source("src/deSelect.R")
source("src/pathEnrich.R")
load("db/pathways.rda")
source("src/ars.R")
source("src/perm.R")
source("src/arsPvalue.R")
load("db/ARSF_345pathway.RData") #ARSFpathway_out: 345x66
load("db/ARSFpvalue_345pathway.RData") #ARSFpvalue: 345x66
source("src/clustPathway.R")
source("src/mdsPathway.R")
source("src/textMine.R")
source("src/mdsModel.R")
source("src/clustModel.R")
source("src/heatmapModel.R")
source("src/heatmapGene.R")
source("src/keggView.R")


source("global/constants.R")
source("global/messages.R")
source("global/help.R")
source("global/database.R")
source("global/helpers.R")
source("global/directoryInput.R")

# Create the directory for database prior to application startup
db <- new("Database", name="studies")

# Include all server modules
dir <- "server"
for (f in list.files(path=dir, pattern="*.R")) {
  source(paste(dir, f, sep="/"))
}

# Include all UI modules
dir <- "ui"
for (f in list.files(path=dir, pattern="*.R")) {
  source(paste(dir, f, sep="/"))
}

# Setting default working sirectory
tryCatch({
  DB.load.working.dir(db)
}, error=function(error){
  DB.set.working.dir(db, paste(getwd(), "data", sep="/"))
})

