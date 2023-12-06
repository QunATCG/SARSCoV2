###
 # @Descripttion: Data Analysis
 # @Author: Qun Li
 # @Email: qun.li@ki.se/liqun95@163.com
 # @Date: 2023-11-23 14:21:57
 # @LastEditTime: 2023-11-23 14:21:58
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
####
# > packageVersion('DESeq2')
# [1] ‘1.34.0’
# > packageVersion('pheatmap')
# [1] ‘1.0.12’
####

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

# exclude duplicated genes


# function definition
qunDEGAnalysisTwoCondition <- function(dat, condition_CT, num_CT, condition_Treat, num_Treat, fc, pthreshold){
  coldata <- data.frame(condition = factor(c(rep(condition_CT, num_CT), rep(condition_Treat, num_Treat)), levels = c(condition_CT, condition_Treat)))
  rownames(coldata) <- colnames(dat)
  dds <- DESeqDataSetFromMatrix(countData = dat, colData = coldata, design= ~condition)
  dds1 <- DESeq(dds, parallel = FALSE)
  res <- results(dds1, contrast = c('condition', condition_Treat, condition_CT))
  res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
  res1 <- res1[order(res1$padj, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
  res1[which(res1$log2FoldChange >= log2(fc) & res1$padj < pthreshold),'sig'] <- 'up'
  res1[which(res1$log2FoldChange <= -log2(fc) & res1$padj < pthreshold),'sig'] <- 'down'
  res1[which(abs(res1$log2FoldChange) <= log2(fc) | res1$padj >= pthreshold),'sig'] <- 'none'
  res1_select <- subset(res1, sig %in% c('up', 'down'))
  resultList <- list(DEG = res1_select, datNormal = dds1, coldata = coldata )
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

filterResult <- function(dat, expData, nameTag, fc, pthreshold, tagItem){
  if(tagItem == "three"){
    res1 <- data.frame(results(dat, name=nameTag), stringsAsFactors = FALSE, check.names = FALSE)
    res1$Ensembl <- rownames(res1)
    res1[which(res1$log2FoldChange >= log2(fc) & res1$padj < pthreshold),'sig'] <- 'up'
    res1[which(res1$log2FoldChange <= -log2(fc) & res1$padj < pthreshold),'sig'] <- 'down'
    res1[which(abs(res1$log2FoldChange) <= log2(fc) | res1$padj >= pthreshold),'sig'] <- 'none'
    res1_select <- subset(res1, sig %in% c('up', 'down'))
    res <- merge(res1_select, expData, by = "Ensembl", all.x = T)
    return(res)
  }
  if(tagItem == "two"){
    res1 <- data.frame(dat, stringsAsFactors = FALSE, check.names = FALSE)
    res1$Ensembl <- rownames(res1)
    res <- merge(res1, expData, by = "Ensembl", all.x = T)
    return(res)
  }
}

# expression filter
# exculde genes with fpkm <= 1 and readcounts <= 5
star_total_over <- keepOverReadCount(keepOverExp(star_total,1),5,"star_total")
star_covid_over <- keepOverReadCount(keepOverExp(star_covid,1),5,"star_covid19")

hisat2_total_over <- keepOverReadCount(keepOverExp(hisat2_total,1),5,"hisat2_total")
hisat2_covid_over <- keepOverReadCount(keepOverExp(hisat2_covid,1),5,"hisat2_covid19")

star_total_over <- star_total
star_covid_over <- star_covid

# PCA
star_total_over_PCA <- star_total_over[,30:38]
rownames(star_total_over_PCA) <- star_total_over$Ensembl
star_total_over_PCA_DEseq2 <- qunDEGAnalysisThreeCondition(round(star_total_over_PCA), "Mock", 3, "NT", 3, "T", 3, 2, 0.01)
plotDispEsts(star_total_over_PCA_DEseq2$DEG, main="Dispersion plot")
rld_star_total_over <- rlog(star_total_over_PCA_DEseq2$DEG)
star_total_over_PCA_dataFrame <- plotPCA(rld_star_total_over, returnData = T)

hisat2_total_over_PCA <- hisat2_total_over[,21:29]
rownames(hisat2_total_over_PCA) <- hisat2_total_over$Ensembl
hisat2_total_over_PCA_DEseq2 <- qunDEGAnalysisThreeCondition(round(hisat2_total_over_PCA), "Mock", 3, "NT", 3, "T", 3, 2, 0.01)
plotDispEsts(hisat2_total_over_PCA_DEseq2$DEG, main="Dispersion plot")
rld_hisat2_total_over <- rlog(hisat2_total_over_PCA_DEseq2$DEG)
hisat2_total_over_PCA_dataFrame <- plotPCA(rld_hisat2_total_over, returnData = T)

#pdf("./Results/0_PCA.pdf", width = 6.29, height = 2.26)
ggplot(data = hisat2_total_over_PCA_dataFrame, mapping = aes(x = PC1, y = PC2, colour = group)) + geom_point(size = 1) +
  scale_color_manual(values = c( Mock = "red", NT = "blue", T = "green" )) +
  xlab("PC1: 90% variance") +
  ylab("PC2: 4% variance") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
#dev.off()

rld_hisat2_total_over_mat <- assay(rld_hisat2_total_over)
rld_hisat2_total_over_cor <- cor(rld_hisat2_total_over_mat)

#pdf("./Results/0_Pheatmap.pdf", width = 6.29, height = 5.68)
pheatmap(rld_hisat2_total_over_cor)
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
resultsNames(star_total_over_PCA_DEseq2$DEG)

#resultsNames(hisat2_total_over_PCA_DEseq2$DEG)

# generate results table for NT_vs_Mock
head(star_total_over_PCA)
DEG_NT_Mock_star <- qunDEGAnalysisTwoCondition(dat = star_total_over_PCA[,1:6],condition_CT = "Mock", num_CT = 3,
                                              condition_Treat = "NT", num_Treat = 3, fc = 2, pthreshold = 0.01)

DEG_NT_Mock_star_result <- filterResult(dat = DEG_NT_Mock_star$DEG, nameTag = "", expData = star_total_over, 
                                        fc = 2, pthreshold = 0.01, tagItem = "two")

DEG_NT_Mock_hisat2 <- qunDEGAnalysisTwoCondition(dat = hisat2_total_over_PCA[,1:6],condition_CT = "Mock", num_CT = 3,
                                               condition_Treat = "NT", num_Treat = 3, fc = 2, pthreshold = 0.01)

DEG_NT_Mock_hisat2_result <- filterResult(dat = DEG_NT_Mock_hisat2$DEG, nameTag = "", expData = hisat2_total_over, 
                                        fc = 2, pthreshold = 0.01, tagItem = "two")

#dat=star_total_over_PCA[,4:9]
#condition_CT="T" 
#num_CT=3
#condition_Treat="NT"
#num_Treat=3
#fc=2
#pthreshold=0.01

DEG_NT_T_star <- qunDEGAnalysisTwoCondition(dat = star_total_over_PCA[,4:9],condition_CT = "NT", num_CT = 3,
                                               condition_Treat = "T", num_Treat = 3, fc = 2, pthreshold = 0.05)

DEG_NT_T_star_result <- filterResult(dat = DEG_NT_T_star$DEG, nameTag = "", expData = star_total_over, 
                                        fc = 2, pthreshold = 0.01, tagItem = "two")


DEG_NT_T_hisat2 <- qunDEGAnalysisTwoCondition(dat = hisat2_total_over_PCA[,4:9],condition_CT = "NT", num_CT = 3,
                                            condition_Treat = "T", num_Treat = 3, fc = 2, pthreshold = 0.05)

DEG_NT_T_hisat2_result <- filterResult(dat = DEG_NT_T_hisat2$DEG, nameTag = "", expData = hisat2_total_over, 
                                     fc = 2, pthreshold = 0.01, tagItem = "two")



testGene <- res_NT_Mock_star_hisat2[res_NT_Mock_star_hisat2$sig == "up",]$Ensembl

#gene.df <- bitr(gene, fromType = "ENTREZID", toType = c("ENSEMBL", "SYMBOL"), OrgDb = org.Hs.eg.db)

ego_ALL <- enrichGO(gene = testGene, 
                    #universe = names(geneList),
                    OrgDb = org.Hs.eg.db,
                    keyType = 'ENSEMBL',
                    ont = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 1,
                    qvalueCutoff = 1,
                    readable = TRUE)

write.table(testGene,"test.gene.txt",quote = F, row.names = F)
# Mock vs NT

