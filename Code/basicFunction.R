# function definition
{
  qunDEGAnalysisTwoCondition <- function(dat, condition_CT, num_CT, condition_Treat, num_Treat){
    coldata <- data.frame(condition = factor(c(rep(condition_CT, num_CT), rep(condition_Treat, num_Treat)), levels = c(condition_CT, condition_Treat)))
    rownames(coldata) <- colnames(dat)
    dds <- DESeqDataSetFromMatrix(countData = dat, colData = coldata, design= ~condition)
    dds1 <- DESeq(dds, parallel = FALSE)
    res <- results(dds1, contrast = c('condition', condition_Treat, condition_CT))
    res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
    res1 <- res1[order(res1$padj, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
    resultList <- list(DEG = res1, datNormal = dds1, coldata = coldata )
    return(resultList)
  }
  
  keepOverExp <- function(dat, cutoff){
    res <- dat[which(rowMeans(dat[,c("Mock_1_fpkm", "Mock_2_fpkm", "Mock_3_fpkm")]) >=cutoff | 
                       rowMeans(dat[,c("T_1_fpkm", "T_2_fpkm", "T_3_fpkm")]) >=cutoff | 
                       rowMeans(dat[,c("NT_1_fpkm", "NT_2_fpkm", "NT_3_fpkm")]) >=cutoff),]
  }
  
  keepOverReadCount <- function(dat, cutoff){
    res <- dat[which(rowSums(dat[,c("Mock_1_readcount", "Mock_2_readcount", "Mock_3_readcount")] >= cutoff ) >= 3 |
                         rowSums(dat[,c("NT_1_readcount", "NT_2_readcount", "NT_3_readcount")] >= cutoff ) >= 3 |
                         rowSums(dat[,c("T_1_readcount", "T_2_readcount", "T_3_readcount")] >= cutoff ) >= 3),]
    return(res)
  }
  
  markResult <- function(dat, fc, pthreshold, tagItem){
    if(tagItem == "treatVScontrol"){
      res1 <- dat[order(dat$pvalue, dat$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
      res1[which(res1$log2FoldChange >= log2(fc) & res1$pvalue < pthreshold),'sig'] <- 'up'
      res1[which(res1$log2FoldChange <= -log2(fc) & res1$pvalue < pthreshold),'sig'] <- 'down'
      res1[which(abs(res1$log2FoldChange) <= log2(fc) | res1$pvalue >= pthreshold),'sig'] <- 'none'
      return(res1)
    }
    if(tagItem == "controlVStreat"){
      res1 <- dat[order(dat$pvalue, dat$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
      res1[which(res1$log2FoldChange >= log2(fc) & res1$pvalue < pthreshold),'sig'] <- 'down'
      res1[which(res1$log2FoldChange <= -log2(fc) & res1$pvalue < pthreshold),'sig'] <- 'up'
      res1[which(abs(res1$log2FoldChange) <= log2(fc) | res1$pvalue >= pthreshold),'sig'] <- 'none'
      return(res1)
    }
  }
  
  qunGO <- function(dat, tagItem){
    dataForGO <- dat
    GeneForGO <- dataForGO[dataForGO$sig == tagItem,]$Ensembl
    ego_ALL <- enrichGO(gene = GeneForGO, 
                        universe = dataForGO$Ensembl,
                        OrgDb = org.Hs.eg.db,
                        keyType = 'ENSEMBL',
                        ont = "ALL",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05,
                        readable = TRUE)
    res <- as.data.frame(ego_ALL)
    return(res)
  }
}