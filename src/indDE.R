indDE <- function(data, group, data.type, case.label, ctrl.label){
  ## data is a numeric matrix, group is a factor
  
  check.compatibility(data, group, case.label, ctrl.label)

  if(data.type=="microarray") {
    #implement limma
    group <- relevel(group,ref=ctrl.label)
    design <-model.matrix(~group)
    fit <-lmFit(data, design)
    ebFit<-eBayes(fit)
    out.table <- topTable(ebFit,coef=2, number=Inf, sort.by='none')
    log2FC <- out.table$logFC
    lfcSE <- sqrt(ebFit$s2.post) * fit$stdev.unscaled[,2]
    p <- as.numeric(out.table$P.Value)
    q <- p.adjust(p,method="BH")
    summary <- data.frame(logFC=log2FC,lfcSE = lfcSE,
                          pvalue = p, qvalue= q)
    rownames(summary) <- rownames(data)
  } 
  
  if(data.type=="RNAseq") {
    #implement DESeq2
    group <- relevel(group,ref=ctrl.level)
    design <-model.matrix(~ group)  # design matrix
    colData <- data.frame(group=group)
    colnames(group) <- colnames(design)[-1] 
    ddsMat <- DESeqDataSetFromMatrix(countData = data,
                                     colData = colData,
                                     design = as.formula(
                            paste(" ~ ",paste(colnames(colData), collapse=" + ") 
                                  ) )  )
    ddsMat <- DESeq(ddsMat)
    res <- results(ddsMat,contrast=c(colnames(colData)[1],levels(group)[2],
                                     levels(group)[1]) )
    log2FC <- as.numeric(res$log2FoldChange)
    lfcSE <- as.numeric(res$lfcSE)
    p <- as.numeric(res$pvalue)
    q <- p.adjust(p,method="BH")
    summary <- data.frame(logFC=log2FC,lfcSE = lfcSE,
                          pvalue = p, qvalue= q)
    rownames(summary) <- rownames(data)
  } 
  return(summary)
}