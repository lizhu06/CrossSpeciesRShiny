read.rawData <- function(data.file, sep=",", quote='"', header=T) {
  #input either csv or txt file, gene on rows and sample on columns
    data <- as.matrix(read.csv(data.file, sep=",", quote=quote, 
                               header=header, row.names=1))
    check.rawData(data)  # check the expression data 
    return(data)
}

read.groupData <- function(group.file, sep=",", quote='"', header=T) {
  #input either csv or txt file, samples on rows, one column containing the class label
  group <- as.factor(read.csv(group.file, sep=",", quote=quote, 
                              header=header, row.names=1)[,1])
  check.groupData(group)  # check the group data
  return(group)
}
  
read.pData <-  function(p.file, sep=",",quote='"', header=T){
  #input either csv or txt file, gene on rows, two columns: 1st: 2-sided p-value, 2nd: logFC (or equivalent effect size)
  pData <- as.matrix(read.csv(p.file, sep=",", quote=quote, 
                              header=header, row.names=1))
  check.pData(pData)  # check the pvalue data 
  return(pData)
} 

