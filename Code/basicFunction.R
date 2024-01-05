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
    res <- dat[which(rowSums(dat[,c("Mock_1_rawcount", "Mock_2_rawcount", "Mock_3_rawcount")] >= cutoff ) >= 3 |
                         rowSums(dat[,c("NT_1_rawcount", "NT_2_rawcount", "NT_3_rawcount")] >= cutoff ) >= 3 |
                         rowSums(dat[,c("T_1_rawcount", "T_2_rawcount", "T_3_rawcount")] >= cutoff ) >= 3),]
    return(res)
  }
  
  markResult <- function(dat, fc, pthreshold, tagItem){
    if(tagItem == "treatVScontrol"){
      res1 <- dat[order(dat$pvalue, dat$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
      res1[which(res1$log2FoldChange >= log2(fc) & res1$pvalue < pthreshold),'sig'] <- 'up'
      res1[which(res1$log2FoldChange <= -log2(fc) & res1$pvalue < pthreshold),'sig'] <- 'down'
      res1[which(abs(res1$log2FoldChange) <= log2(fc) | res1$pvalue >= pthreshold),'sig'] <- 'none'
      # strict
      res1[which(res1$log2FoldChange >= log2(2) & res1$padj < 0.05),'sig_strict'] <- 'up'
      res1[which(res1$log2FoldChange <= -log2(2) & res1$padj < 0.05),'sig_strict'] <- 'down'
      res1[which(abs(res1$log2FoldChange) <= log2(2) | res1$padj >= 0.05),'sig_strict'] <- 'none'
      return(res1)
    }
    if(tagItem == "controlVStreat"){
      res1 <- dat[order(dat$pvalue, dat$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
      res1[which(res1$log2FoldChange >= log2(fc) & res1$pvalue < pthreshold),'sig'] <- 'down'
      res1[which(res1$log2FoldChange <= -log2(fc) & res1$pvalue < pthreshold),'sig'] <- 'up'
      res1[which(abs(res1$log2FoldChange) <= log2(fc) | res1$pvalue >= pthreshold),'sig'] <- 'none'
      # strict
      res1[which(res1$log2FoldChange >= log2(2) & res1$padj < 0.05),'sig_strict'] <- 'down'
      res1[which(res1$log2FoldChange <= -log2(2) & res1$padj < 0.05),'sig_strict'] <- 'up'
      res1[which(abs(res1$log2FoldChange) <= log2(2) | res1$padj >= 0.05),'sig_strict'] <- 'none'
      return(res1)
    }
  }
  
  qunGO <- function(dat, tagItem, isStrict){
    dataForGO <- dat
    if(isStrict == "T"){
      GeneForGO <- dataForGO[dataForGO$sig_strict == tagItem,]$Ensembl
    }
    if(isStrict == "N"){
      GeneForGO <- dataForGO[dataForGO$sig == tagItem,]$Ensembl
    }
    ego_ALL <- enrichGO(gene = GeneForGO, 
                        universe = dataForGO$Ensembl,
                        OrgDb = org.Hs.eg.db,
                        keyType = 'ENSEMBL',
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.1,
                        readable = TRUE)
    res1 <- as.data.frame(ego_ALL)
    res <- list(ego = ego_ALL, result = res1)
    return(res)
  }
  
  qunplotValcano <- function(dat, tagItem, rangeNum){
    #dat <- na.omit(dat[,-7])
    ggplot(dat,aes(log2FoldChange,-log10(pvalue),color = sig))+ 
      geom_point()+
      scale_x_continuous(limits = c(-rangeNum, rangeNum)) +
      scale_color_manual(values = c(down = "Dark Green", up = "Dark Red", none = "grey")) +
      labs(x= expression(Log[2]*" Fold Change"), y = expression(-Log[10]*" (pvalue)"), 
           title = paste(tagItem, "up: ", table(dat$sig)[3], "Down: ", table(dat$sig)[1], sep = " ")) +
      theme_classic()
  }
  
  qunplotGO <- function(dat){
    ggplot(dat, aes(-log10(pvalue), reorder(Description, -log10(pvalue), decreasing = T))) +
      geom_bar(stat = "identity") + 
      xlab("-log10(pvalue)") + ylab("GO term") + 
      theme_bw() + theme_classic()
  }
  
  qunNetworkGeneration <- function(dat){
    net_dataframe <- data.frame(matrix(nrow = cumsum(1:nrow(dat))[nrow(dat)], ncol = 3))
    colnames(net_dataframe) <- c("GOID1", "GOID2", "CommonNum")
    countN <- 0
    for (i in 1:nrow(dat)){
      nodei_ID <- dat[i,]$ID
      nodei_Gene <- strsplit(dat[i,]$geneID, "/")[[1]]
      for (j in i:nrow(dat)){
        nodej_ID <- dat[j,]$ID
        nodej_Gene <- strsplit(dat[j,]$geneID, "/")[[1]]
        num_inter <- length(intersect(nodei_Gene, nodej_Gene))
        countN <- countN + 1
        net_dataframe[countN,1] <- nodei_ID
        net_dataframe[countN,2] <- nodej_ID
        net_dataframe[countN,3] <- num_inter
      }
    }
    net_dataframe <- net_dataframe[net_dataframe[,1] != net_dataframe[,2],]
    net_dataframe <- net_dataframe[net_dataframe[,3] != 0,]
    node_dataframe <- dat[,c("ID", "Description","pvalue", "Count", "geneID")]
    res <- list(net = net_dataframe, node = node_dataframe)
    return(res)
  }
  
  qunplotbox <- function(dat_net, dat_exp, itemType, myComparision, myColors,isThisStudy){
    dat_net_select <- dat_net[dat_net$Type == itemType,]$geneID
    if (length(dat_net_select) == 1){
      genes_select <- strsplit(dat_net_select[1], "/")[[1]]
    }else{
      genes_select_tmp <- strsplit(dat_net_select[1], "/")[[1]]
      for (i in 2:length(dat_net_select)){
        genes_select_i <- strsplit(dat_net_select[i], "/")[[1]]
        genes_select_tmp <- union(genes_select_tmp, genes_select_i)
      }
      genes_select <- genes_select_tmp
    }
    if (isThisStudy == "T"){
      genes_select_exp <- dat_exp[dat_exp$Gene %in% genes_select,][,c("Mock_mean_scale", "NT_mean_scale", "T_mean_scale")]
      colnames(genes_select_exp) <- c("Mock", "NT", "T")
      genes_select_exp_stack <- stack(genes_select_exp)
    }
    if (isThisStudy == "F"){
      genes_select_exp <- dat_exp[dat_exp$Gene %in% genes_select,][,c("NC_Mean_scale", "Case_Mean_scale")]
      colnames(genes_select_exp) <- c("NC", "Case")
      genes_select_exp_stack <- stack(genes_select_exp)
    }
    # myComparision <- list(c("Mock_mean_scale", "NT_mean_scale"))
    ggplot(genes_select_exp_stack, aes(x=ind, y=values, fill = ind)) + 
      geom_boxplot(outlier.shape = NA) + 
      scale_fill_manual(values = myColors) +
      #stat_summary(fun.data = mean_sdl,fun.args = list(mult=1), geom = "errorbar", width = 0.2)+
      #stat_compare_means(comparisons = myComparision, method = "wilcox.test", label = "p.forma") +
      labs(x = "",y = "Average Expression level (Z scaled)") +
      theme_classic() + theme(legend.position="none")
      #theme(axis.text.x = element_blank())
  }
  
  qunplotboxgenes <- function(gene, dat_exp, myComparision, myColors,isThisStudy){
    if (isThisStudy == "T"){
      gene_select_exp <- dat_exp[dat_exp$Gene %in% gene,][,c("Mock_1_fpkm", "Mock_2_fpkm", "Mock_3_fpkm",
                                                                          "NT_1_fpkm", "NT_2_fpkm", "NT_3_fpkm", "T_1_fpkm", 
                                                                          "T_2_fpkm", "T_3_fpkm")]
      colnames(gene_select_exp) <- c(rep("Mock", 3), rep("NT", 3), rep("T", 3))
      gene_select_exp_stack <- stack(gene_select_exp)
      gene_select_exp_stack$ind <- c(rep("Mock", 3), rep("NT", 3), rep("T", 3))
    }
    if (isThisStudy == "F"){
      tmp_colnames <- c("Ensembl", "Gene")
      colnames_gene_select_exp <- colnames(dat_exp)[grep("Mean",colnames(dat_exp), invert = TRUE)]
      colnames_gene_select_exp <- setdiff(colnames_gene_select_exp, tmp_colnames)
      colnames_gene_select_exp_case <- colnames_gene_select_exp[grep("Case", colnames_gene_select_exp)]
      colnames_gene_select_exp_nc <- colnames_gene_select_exp[grep("NC", colnames_gene_select_exp)]
      gene_select_exp <- dat_exp[dat_exp$Gene %in% gene,][,c(colnames_gene_select_exp_nc, colnames_gene_select_exp_case)]
      
      colnames(gene_select_exp) <- c(rep("NC", length(colnames_gene_select_exp_nc)), rep("Case", length(colnames_gene_select_exp_case)))
      gene_select_exp_stack <- stack(gene_select_exp)
      #print(gene_select_exp_stack)
      gene_select_exp_stack$ind <- c(rep("NC", length(colnames_gene_select_exp_nc)), rep("Case", length(colnames_gene_select_exp_case)))
      #print(gene_select_exp_stack)
      gene_select_exp_stack$ind <- factor(x = gene_select_exp_stack$ind, levels = c("NC", "Case"))
    }
    # myComparision <- list(c("Mock_mean_scale", "NT_mean_scale"))
    ggplot(gene_select_exp_stack, aes(x=ind, y=values, fill = ind)) + 
      geom_boxplot(outlier.shape = NA) + 
      scale_fill_manual(values = myColors) +
      #stat_summary(fun.data = mean_sdl,fun.args = list(mult=1), geom = "errorbar", width = 0.2)+
      #stat_compare_means(comparisons = myComparision, method = "wilcox.test", label = "p.forma") +
      labs(x = "",y = "Average Expression level (Z scaled)") +
      theme_classic()
  }
  
  qunplotmultiupbox <- function(dat_net, exp1, exp2, exp3, exp4, itemType){
    print("this study")
    p_box_thisStudy <- qunplotbox(dat_net, dat_exp = exp1, itemType, 
                                                    myComparision = list(c("Mock", "NT]"), c("Mock", "T"), c("NT", "T")), 
                                                    myColors = c("grey", "Dark Red", "Brown"), isThisStudy = "T")
    
    print("GSE152418")
    p_box_GSE152418 <- qunplotbox(dat_net, dat_exp = exp2, itemType = itemType,
                                         myComparision = list(c("NC", "Case")),
                                         myColors = c("grey", "Dark Red"), isThisStudy = "F")
    
    print("GSE157103")
    p_box_GSE157103 <- qunplotbox(dat_net, dat_exp = exp3, itemType = itemType,
                                         myComparision = list(c("NC", "Case")),
                                         myColors = c("grey", "Dark Red"), isThisStudy = "F")
    
    print("GSE150316")
    p_box_GSE150316 <- qunplotbox(dat_net, dat_exp = exp4, itemType = itemType,
                                  myComparision = list(c("NC", "Case")),
                                  myColors = c("grey", "Dark Red"), isThisStudy = "F")
    
    res <- list(p_thisStudy = p_box_thisStudy, p_GSE152418 = p_box_GSE152418,
                p_GSE157103 = p_box_GSE157103, p_GSE150316 = p_box_GSE150316)
    return(res)
  }
  
  qunplotmultidownbox <- function(dat_net, exp1, exp2, exp3, exp4, itemType){
    print("this study")
    p_box_thisStudy <- qunplotbox(dat_net = dat_net, dat_exp = exp1, itemType = itemType, 
                                  myComparision = list(c("Mock", "NT"), c("Mock", "T"), c("NT", "T")), 
                                  myColors = c("black","Dark Green", "Olive Drab"), isThisStudy = "T")
    print("GSE152418")
    p_box_GSE152418 <- qunplotbox(dat_net, dat_exp = exp2, itemType = itemType,
                                  myComparision = list(c("NC", "Case")),
                                  myColors = c("black","Dark Green"), isThisStudy = "F")
    print("GSE157103")
    p_box_GSE157103 <- qunplotbox(dat_net, dat_exp = exp3, itemType = itemType,
                                  myComparision = list(c("NC", "Case")),
                                  myColors = c("black","Dark Green"), isThisStudy = "F")
    print("GSE150316")
    p_box_GSE150316 <- qunplotbox(dat_net, dat_exp = exp4, itemType = itemType,
                                  myComparision = list(c("NC", "Case")),
                                  myColors = c("black","Dark Green"), isThisStudy = "F")
    
    res <- list(p_thisStudy = p_box_thisStudy, p_GSE152418 = p_box_GSE152418,
                p_GSE157103 = p_box_GSE157103, p_GSE150316 = p_box_GSE150316)
    return(res)
  }
  
  qunloadDepth <- function(filePath, colnamesTag){
    depth_data <- read.table(filePath, header = FALSE, sep = "\t")
    colnames(depth_data) <- c("Chr", "Star", "End", paste("Coverage", colnamesTag, sep = "_"))
    depth_data$ID <- paste(depth_data$Star, depth_data$End, sep = "_")
    return(depth_data)
  }
}
