###
 # @Descripttion: Data Analysis
 # @Author: Qun Li
 # @Email: qun.li@ki.se/liqun95@163.com
 # @Date: 2023-11-23 14:21:57
### 

# DEG analysis was performed using R package DESeq2
# ploting was performed using R package pheatmap

# set working env
rm(list = ls())
# set your working path 
setwd("/Users/liqun/Desktop/Projects/Covid19/Data/Code/SARSCoV2/")
# check your working path
dirNow <- getwd()

# load libraries
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

# load data
star_total <- read.table("./Data/ExpRNAseq/STAR/Star_exp_total.txt", header = T,
                         sep = "\t", stringsAsFactors = FALSE)
star_covid <- read.table("./Data/ExpRNAseq/STAR/Star_exp_Covid19.txt", header = T,
                         sep = "\t", stringsAsFactors = FALSE)

hisat2_total <- read.table("./Data/ExpRNAseq/Hisat2/Hisat2_exp_total.txt", header = T,
                           sep = "\t", stringsAsFactors = FALSE)
hisat2_covid <- read.table("./Data/ExpRNAseq/Hisat2/Hisat2_exp_Covid19.txt", header = T,
                           sep = "\t", stringsAsFactors = FALSE)

Covid19_geneName <- c("ORF1ab","ORF1ab","S", "ORF3a", "E", "M", "ORF6", "ORF7a","ORF7b",
                      "ORF8", "N", "ORF10")
star_covid$Gene <- Covid19_geneName
hisat2_covid$Gene <- Covid19_geneName

# function definition
{
  qunDEGAnalysisTwoCondition <- function(dat, condition_CT, num_CT, condition_Treat, num_Treat, fc, pthreshold){
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
  
  qunDEGAnalysisThreeCondition <- function(dat, condition_CT, num_CT, condition_Treat1, num_Treat1, condition_Treat2, num_Treat2, fc, pthreshold){
    coldata <- data.frame(condition = factor(c(rep(condition_CT, num_CT), rep(condition_Treat1, num_Treat1), rep(condition_Treat2, num_Treat2)), levels = c(condition_CT, condition_Treat1, condition_Treat2)))
    rownames(coldata) <- colnames(dat)
    dds <- DESeqDataSetFromMatrix(countData = dat, colData = coldata, design= ~condition)
    dds1 <- DESeq(dds, parallel = FALSE)
    resultList <- list(DEG =  dds1, coldata = coldata )
    return(resultList)
  }
  
  keepOverExp <- function(dat, cutoff){
    res <- dat[which(rowMeans(dat[,c("Mock_1_FPKM", "Mock_2_FPKM", "Mock_3_FPKM")]) >=cutoff | 
                       rowMeans(dat[,c("T_1_FPKM", "T_2_FPKM", "T_3_FPKM")]) >=cutoff | 
                       rowMeans(dat[,c("NT_1_FPKM", "NT_2_FPKM", "NT_3_FPKM")]) >=cutoff),]
  }
  
  keepOverReadCount <- function(dat, cutoff, tagItem){
    if(tagItem == "star_total"){
      res <- dat[which(rowSums(dat[,c("Mock_1_ReadCount", "Mock_2_ReadCount", "Mock_3_ReadCount")] >= cutoff ) >= 3 |
                         rowSums(dat[,c("NT_1_ReadCount", "NT_2_ReadCount", "NT_3_ReadCount")] >= cutoff ) >= 3 |
                         rowSums(dat[,c("T_1_ReadCount", "T_2_ReadCount", "T_3_ReadCount")] >= cutoff ) >= 3),]
    }
    if(tagItem == "star_covid19"){
      res <- dat[which(rowSums(dat[,c("Mock_1_ExpectedCount", "Mock_2_ExpectedCount", "Mock_3_ExpectedCount")] >= cutoff ) >= 3 |
                         rowSums(dat[,c("NT_1_ExpectedCount", "NT_2_ExpectedCount", "NT_3_ExpectedCount")] >= cutoff ) >= 3 |
                         rowSums(dat[,c("T_1_ExpectedCount", "T_2_ExpectedCount", "T_3_ExpectedCount")] >= cutoff ) >= 3),]
    }
    if(tagItem == "hisat2_total"){
      res <- dat[which(rowSums(dat[,c("Mock_1_ReadCount", "Mock_2_ReadCount", "Mock_3_ReadCount")] >= cutoff ) >= 3 |
                         rowSums(dat[,c("NT_1_ReadCount", "NT_2_ReadCount", "NT_3_ReadCount")] >= cutoff ) >= 3 |
                         rowSums(dat[,c("T_1_ReadCount", "T_2_ReadCount", "T_3_ReadCount")] >= cutoff ) >= 3),]
    }
    if(tagItem == "hisat2_covid19"){
      res <- dat[which(rowSums(dat[,c("Mock_1_ReadCount", "Mock_2_ReadCount", "Mock_3_ReadCount")] >= cutoff ) >= 3 |
                         rowSums(dat[,c("NT_1_ReadCount", "NT_2_ReadCount", "NT_3_ReadCount")] >= cutoff ) >= 3 |
                         rowSums(dat[,c("T_1_ReadCount", "T_2_ReadCount", "T_3_ReadCount")] >= cutoff ) >= 3),]
    }
    return(res)
  }
  
  mergeResult <- function(dat, expData, nameTag, fc, pthreshold, tagItem){
    if(tagItem == "three"){
      res1 <- data.frame(results(dat, name=nameTag), stringsAsFactors = FALSE, check.names = FALSE)
      res1$Ensembl <- rownames(res1)
      res <- merge(res1, expData, by = "Ensembl", all.x = T)
      return(res)
    }
    if(tagItem == "two"){
      res1 <- data.frame(dat, stringsAsFactors = FALSE, check.names = FALSE)
      res1$Ensembl <- rownames(res1)
      res <- merge(res1, expData, by = "Ensembl", all.x = T)
      return(res)
    }
  }
  
  markResult <- function(dat, fc, pthreshold, tagItem){
    if(tagItem == "treatVScontrol"){
      res1 <- dat[order(dat$padj, dat$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
      res1[which(res1$log2FoldChange >= log2(fc) & res1$padj < pthreshold),'sig'] <- 'up'
      res1[which(res1$log2FoldChange <= -log2(fc) & res1$padj < pthreshold),'sig'] <- 'down'
      res1[which(abs(res1$log2FoldChange) <= log2(fc) | res1$padj >= pthreshold),'sig'] <- 'none'
      return(res1)
    }
    if(tagItem == "controlVStreat"){
      res1 <- dat[order(dat$padj, dat$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
      res1[which(res1$log2FoldChange >= log2(fc) & res1$padj < pthreshold),'sig'] <- 'down'
      res1[which(res1$log2FoldChange <= -log2(fc) & res1$padj < pthreshold),'sig'] <- 'up'
      res1[which(abs(res1$log2FoldChange) <= log2(fc) | res1$padj >= pthreshold),'sig'] <- 'none'
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


# expression filter
# exculde genes with fpkm <= 1 and readcounts <= 5
# star_total_over <- keepOverReadCount(keepOverExp(star_total,1),5,"star_total")
# star_covid_over <- keepOverReadCount(keepOverExp(star_covid,1),5,"star_covid19")

# hisat2_total_over <- keepOverReadCount(keepOverExp(hisat2_total,1),5,"hisat2_total")
# hisat2_covid_over <- keepOverReadCount(keepOverExp(hisat2_covid,1),5,"hisat2_covid19")

star_total_over <- star_total
hisat2_total_over <- hisat2_total

# PCA
star_total_over_PCA <- star_total_over[,30:38]
rownames(star_total_over_PCA) <- star_total_over$Ensembl
# star_total_over_PCA_DEseq2 <- qunDEGAnalysisThreeCondition(round(star_total_over_PCA), "Mock", 3, "NT", 3, "T", 3, 2, 0.01)
#plotDispEsts(star_total_over_PCA_DEseq2$DEG, main="Dispersion plot")
# rld_star_total_over <- rlog(star_total_over_PCA_DEseq2$DEG)
# star_total_over_PCA_dataFrame <- plotPCA(rld_star_total_over, returnData = T)

hisat2_total_over_PCA <- hisat2_total_over[,21:29]
rownames(hisat2_total_over_PCA) <- hisat2_total_over$Ensembl
#hisat2_total_over_PCA_DEseq2 <- qunDEGAnalysisThreeCondition(round(hisat2_total_over_PCA), "Mock", 3, "NT", 3, "T", 3, 2, 0.01)
#plotDispEsts(hisat2_total_over_PCA_DEseq2$DEG, main="Dispersion plot")
#rld_hisat2_total_over <- rlog(hisat2_total_over_PCA_DEseq2$DEG)
#hisat2_total_over_PCA_dataFrame <- plotPCA(rld_hisat2_total_over, returnData = T)

#pdf("./Results/0_PCA.pdf", width = 6.29, height = 2.26)
#ggplot(data = hisat2_total_over_PCA_dataFrame, mapping = aes(x = PC1, y = PC2, colour = group)) + geom_point(size = 1) +
#  scale_color_manual(values = c( Mock = "red", NT = "blue", T = "green" )) +
#  xlab("PC1: 90% variance") +
#  ylab("PC2: 4% variance") +
#  theme_bw() +
#  theme(panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank())
#dev.off()

#rld_hisat2_total_over_mat <- assay(rld_hisat2_total_over)
#rld_hisat2_total_over_cor <- cor(rld_hisat2_total_over_mat)

#pdf("./Results/0_Pheatmap.pdf", width = 6.29, height = 5.68)
#pheatmap(rld_hisat2_total_over_cor)
#dev.off()

#vst_star_total_over1 <- vst(star_total_over1_PCA_DEseq2$DEG)
#plotPCA(vst_star_total_over1)

#hisat2_total_over_PCA <- hisat2_total_over[,21:29]
#rownames(hisat2_total_over_PCA) <- hisat2_total_over$Ensembl
#hisat2_total_over_PCA_DEseq2 <- qunDEGAnalysisThreeCondition(round(hisat2_total_over_PCA), "Mock", 3, "NT", 3, "T", 3, 2, 0.01)
#plotDispEsts(hisat2_total_over1_PCA_DEseq2$DEG, main="Dispersion plot")
#rld_hisat2_total_over1 <- rlog(hisat2_total_over1_PCA_DEseq2$DEG)
#plotPCA(rld_hisat2_total_over1)
#rld_hisat2_total_over1_mat <- assay(rld_hisat2_total_over1)
#rld_hisat2_total_over1_cor <- cor(rld_hisat2_total_over1_mat) 
#pheatmap(rld_hisat2_total_over1_cor)
#vst_hisat2_total_over1 <- vst(hisat2_total_over1_PCA_DEseq2$DEG)
#plotPCA(vst_hisat2_total_over1)

# DEG analysis
#resultsNames(star_total_over_PCA_DEseq2$DEG)

#resultsNames(hisat2_total_over_PCA_DEseq2$DEG)

# generate results table for NT_vs_Mock
# head(star_total_over_PCA)
DEG_NT_Mock_star <- qunDEGAnalysisTwoCondition(dat = star_total_over_PCA[,1:6],condition_CT = "Mock", num_CT = 3,
                                              condition_Treat = "NT", num_Treat = 3, fc = 2, pthreshold = 0.01)

DEG_NT_Mock_star_result <- markResult(dat = mergeResult(dat = DEG_NT_Mock_star$DEG, nameTag = "", expData = star_total_over, 
                                        fc = 2, pthreshold = 0.01, tagItem = "two"), fc = 1.5, pthreshold = 0.05,tagItem = "treatVScontrol")

DEG_NT_Mock_hisat2 <- qunDEGAnalysisTwoCondition(dat = hisat2_total_over_PCA[,1:6],condition_CT = "Mock", num_CT = 3,
                                               condition_Treat = "NT", num_Treat = 3, fc = 2, pthreshold = 0.01)

DEG_NT_Mock_hisat2_result <- markResult(dat = mergeResult(dat = DEG_NT_Mock_hisat2$DEG, nameTag = "", expData = hisat2_total_over, 
                                        fc = 2, pthreshold = 0.01, tagItem = "two"), fc = 1.5, pthreshold = 0.05,tagItem = "treatVScontrol")



# generate results table for NT_vs_N
DEG_NT_T_star <- qunDEGAnalysisTwoCondition(dat = star_total_over_PCA[,4:9],condition_CT = "NT", num_CT = 3,
                                               condition_Treat = "T", num_Treat = 3, fc = 2, pthreshold = 0.01)

DEG_NT_T_star_result <- markResult(dat = mergeResult(dat = DEG_NT_T_star$DEG, nameTag = "", expData = star_total_over, 
                                        fc = 2, pthreshold = 0.01, tagItem = "two"), fc = 1.5, pthreshold = 0.05,tagItem = "treatVScontrol")


DEG_NT_T_hisat2 <- qunDEGAnalysisTwoCondition(dat = hisat2_total_over_PCA[,4:9],condition_CT = "NT", num_CT = 3,
                                            condition_Treat = "T", num_Treat = 3, fc = 2, pthreshold = 0.01)

DEG_NT_T_hisat2_result <- markResult(dat = mergeResult(dat = DEG_NT_T_hisat2$DEG, nameTag = "", expData = hisat2_total_over, 
                                     fc = 2, pthreshold = 0.01, tagItem = "two"), fc = 1.5, pthreshold = 0.05,tagItem = "treatVScontrol")


# generate results table for T_vs_Mock
DEG_T_Mock_star <- qunDEGAnalysisTwoCondition(dat = star_total_over_PCA[,c(1,2,3,7,8,9)], condition_CT = "Mock", num_CT = 3,
                                               condition_Treat = "T", num_Treat = 3, fc = 2, pthreshold = 0.01)

DEG_T_Mock_star_result <- markResult(dat = mergeResult(dat = DEG_T_Mock_star$DEG, nameTag = "", expData = star_total_over, 
                                        fc = 2, pthreshold = 0.01, tagItem = "two"), fc = 1.5, pthreshold = 0.05,tagItem = "treatVScontrol")

DEG_T_Mock_hisat2 <- qunDEGAnalysisTwoCondition(dat = hisat2_total_over_PCA[,c(1,2,3,7,8,9)],condition_CT = "Mock", num_CT = 3,
                                                 condition_Treat = "T", num_Treat = 3, fc = 2, pthreshold = 0.01)

DEG_T_Mock_hisat2_result <- markResult(dat = mergeResult(dat = DEG_T_Mock_hisat2$DEG, nameTag = "", expData = hisat2_total_over, 
                                          fc = 2, pthreshold = 0.01, tagItem = "two"), fc = 1.5, pthreshold = 0.05,tagItem = "treatVScontrol")



# save DEG result
write.csv2(DEG_NT_Mock_hisat2_result, "./Results/Table/DEG/Hisat2/DEG_NT_Mock_hisat2_result.csv", 
            row.names = FALSE,quote = FALSE)
write.csv2(DEG_NT_Mock_star_result, "./Results/Table/DEG/Star/DEG_NT_Mock_star_result.csv", 
           row.names = FALSE,quote = FALSE)

write.csv2(DEG_T_Mock_hisat2_result, "./Results/Table/DEG/Hisat2/DEG_T_Mock_hisat2_result.csv", 
           row.names = FALSE,quote = FALSE)
write.csv2(DEG_T_Mock_star_result, "./Results/Table/DEG/Star/DEG_T_Mock_star_result.csv", 
           row.names = FALSE,quote = FALSE)

write.csv2(DEG_NT_T_hisat2_result, "./Results/Table/DEG/Hisat2/DEG_NT_T_hisat2_result.csv", 
           row.names = FALSE,quote = FALSE)
write.csv2(DEG_NT_T_star_result, "./Results/Table/DEG/Star/DEG_NT_T_star_result.csv", 
           row.names = FALSE,quote = FALSE)


###############
# DEG_NT_Mock_hisat2_result
# DEG_NT_Mock_star_result
# DEG_T_Mock_hisat2_result
# DEG_T_Mock_star_result
# DEG_NT_T_hisat2_result
# DEG_NT_T_star_result
###############

# GO analysis
#gene.df <- bitr(gene, fromType = "ENTREZID", toType = c("ENSEMBL", "SYMBOL"), OrgDb = org.Hs.eg.db)
# hisat2
DEG_NT_Mock_hisat2_GO_Up <- qunGO(keepOverExp(keepOverReadCount(DEG_NT_Mock_hisat2_result, 5, "hisat2_total"),1), "up")
DEG_NT_Mock_hisat2_GO_Down <- qunGO(keepOverExp(keepOverReadCount(DEG_NT_Mock_hisat2_result, 5, "hisat2_total"),1), "down")

DEG_NT_T_hisat2_GO_Up <- qunGO(keepOverExp(keepOverReadCount(DEG_NT_T_hisat2_result, 5, "hisat2_total"),1), "up")
DEG_NT_T_hisat2_GO_Down <- qunGO(keepOverExp(keepOverReadCount(DEG_NT_T_hisat2_result, 5, "hisat2_total"),1), "down")

DEG_T_Mock_hisat2_GO_Up <- qunGO(keepOverExp(keepOverReadCount(DEG_T_Mock_hisat2_result, 5, "hisat2_total"),1), "up")
DEG_T_Mock_hisat2_GO_Down <- qunGO(keepOverExp(keepOverReadCount(DEG_T_Mock_hisat2_result, 5, "hisat2_total"),1), "down")

# star
DEG_NT_Mock_star_GO_Up <- qunGO(keepOverExp(keepOverReadCount(DEG_NT_Mock_star_result, 5, "star_total"),1), "up")
DEG_NT_Mock_star_GO_Down <- qunGO(keepOverExp(keepOverReadCount(DEG_NT_Mock_star_result, 5, "star_total"),1), "down")

DEG_NT_T_star_GO_Up <- qunGO(keepOverExp(keepOverReadCount(DEG_NT_T_star_result, 5, "star_total"),1), "up")
DEG_NT_T_star_GO_Down <- qunGO(keepOverExp(keepOverReadCount(DEG_NT_T_star_result, 5, "star_total"),1), "down")

DEG_T_Mock_star_GO_Up <- qunGO(keepOverExp(keepOverReadCount(DEG_T_Mock_star_result, 5, "star_total"),1), "up")
DEG_T_Mock_star_GO_Down <- qunGO(keepOverExp(keepOverReadCount(DEG_T_Mock_star_result, 5, "star_total"),1), "down")

#### save DEG results
write.csv2(DEG_NT_Mock_hisat2_GO_Up, "./Results/Table/GO/Hisat2/DEG_NT_Mock_hisat2_GO_Up.csv", row.names = FALSE,quote = FALSE)
write.csv2(DEG_NT_Mock_hisat2_GO_Down, "./Results/Table/GO/Hisat2/DEG_NT_Mock_hisat2_GO_Down.csv", row.names = FALSE,quote = FALSE)
write.csv2(DEG_NT_T_hisat2_GO_Up, "./Results/Table/GO/Hisat2/DEG_NT_T_hisat2_GO_Up.csv", row.names = FALSE,quote = FALSE)
write.csv2(DEG_NT_T_hisat2_GO_Down, "./Results/Table/GO/Hisat2/DEG_NT_T_hisat2_GO_Down.csv", row.names = FALSE,quote = FALSE)
write.csv2(DEG_T_Mock_hisat2_GO_Up, "./Results/Table/GO/Hisat2/DEG_T_Mock_hisat2_GO_Up.csv", row.names = FALSE,quote = FALSE)
write.csv2(DEG_T_Mock_hisat2_GO_Down, "./Results/Table/GO/Hisat2/DEG_T_Mock_hisat2_GO_Down.csv", row.names = FALSE,quote = FALSE)

write.csv2(DEG_NT_Mock_star_GO_Up, "./Results/Table/GO/Star/DEG_NT_Mock_star_GO_Up.csv", row.names = FALSE,quote = FALSE)
write.csv2(DEG_NT_Mock_star_GO_Down, "./Results/Table/GO/Star/DEG_NT_Mock_star_GO_Down.csv", row.names = FALSE,quote = FALSE)
write.csv2(DEG_NT_T_star_GO_Up, "./Results/Table/GO/Star/DEG_NT_T_star_GO_Up.csv", row.names = FALSE,quote = FALSE)
write.csv2(DEG_NT_T_star_GO_Down, "./Results/Table/GO/Star/DEG_NT_T_star_GO_Doan.csv", row.names = FALSE,quote = FALSE)
write.csv2(DEG_T_Mock_star_GO_Up, "./Results/Table/GO/Star/DEG_T_Mock_star_GO_Up.csv", row.names = FALSE,quote = FALSE)
write.csv2(DEG_T_Mock_star_GO_Down, "./Results/Table/GO/Star/DEG_T_Mock_star_GO_Down.csv", row.names = FALSE,quote = FALSE)


# test
#dataForGO <-keepOverReadCount(DEG_NT_Mock_hisat2_result, 5, "hisat2_total")
#GeneForGO <- dataForGO[dataForGO$sig == "up",]$Ensembl
#ego_ALL <- enrichGO(gene = GeneForGO, 
#                    universe = dataForGO$Ensembl,
#                    OrgDb = org.Hs.eg.db,
#                    keyType = 'ENSEMBL',
#                    ont = "ALL",
#                    pAdjustMethod = "BH",
#                   pvalueCutoff = 0.05,
#                   qvalueCutoff = 0.05,
#                    readable = TRUE)
#summary(ego_ALL)


