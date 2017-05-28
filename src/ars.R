ars <- function(dat1,dat2,deIndex1,deIndex2,
                measure=c("youden","Fmeasure","geo.mean")){
  ## same for both global and pathway
  if(measure=="youden"){
     ARS <- ARSY(dat1,dat2,deIndex1)
  } else if(measure=="Fmeasure"){
     ARS <- ARSF(dat1,dat2,deIndex1,deIndex2)
  } else if(measure=="geo.mean"){
     ARS <- ARSG(dat1,dat2,deIndex1)
  }
  return(ARS)
}

