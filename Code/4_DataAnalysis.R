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

keepOver1 <- function(dat){
  res <- dat[which(rowMeans(dat[,c("Mock_1_FPKM", "Mock_2_FPKM", "Mock_3_FPKM")]) >=1 | 
                         rowMeans(dat[,c("T_1_FPKM", "T_2_FPKM", "T_3_FPKM")]) >=1 | 
                         rowMeans(dat[,c("NT_1_FPKM", "NT_2_FPKM", "NT_3_FPKM")]) >=1),]
}
# expression filter
# exculde genes with fpkm <= 1 
star_total_over1 <- keepOver1(star_total)
star_covid_over1 <- keepOver1(star_covid)

hisat2_total_over1 <- keepOver1(hisat2_total)
hisat2_covid_over1 <- keepOver1(hisat2_covid)


# PCA
star_total_over1_PCA <- star_total_over1[,30:38]
rownames(star_total_over1_PCA) <- star_total_over1$Ensembl
star_total_over1_PCA_DEseq2 <- qunDEGAnalysisThreeCondition(round(star_total_over1_PCA), "Mock", 3, "NT", 3, "T", 3, 2, 0.01)
plotDispEsts(star_total_over1_PCA_DEseq2$DEG, main="Dispersion plot")
rld_star_total_over1 <- rlog(star_total_over1_PCA_DEseq2$DEG)
plotPCA(rld_star_total_over1)
rld_star_total_over1_mat <- assay(rld_star_total_over1)
rld_star_total_over1_cor <- cor(rld_star_total_over1_mat) 
pheatmap(rld_star_total_over1_cor)

vst_star_total_over1 <- vst(star_total_over1_PCA_DEseq2$DEG)
plotPCA(vst_star_total_over1)

hisat2_total_over1_PCA <- hisat2_total_over1[,21:29]
rownames(hisat2_total_over1_PCA) <- hisat2_total_over1$Ensembl
hisat2_total_over1_PCA_DEseq2 <- qunDEGAnalysisThreeCondition(round(hisat2_total_over1_PCA), "Mock", 3, "NT", 3, "T", 3, 2, 0.01)
plotDispEsts(hisat2_total_over1_PCA_DEseq2$DEG, main="Dispersion plot")
rld_hisat2_total_over1 <- rlog(hisat2_total_over1_PCA_DEseq2$DEG)
plotPCA(rld_hisat2_total_over1)
rld_hisat2_total_over1_mat <- assay(rld_hisat2_total_over1)
rld_hisat2_total_over1_cor <- cor(rld_hisat2_total_over1_mat) 
pheatmap(rld_hisat2_total_over1_cor)

vst_hisat2_total_over1 <- vst(hisat2_total_over1_PCA_DEseq2$DEG)
plotPCA(vst_hisat2_total_over1)



# 
resultsNames(mm$DEG)

# generate results table for NT_vs_Mock
filterResult <- function(dat, nameTag, fc, pthreshold){
  res1 <- data.frame(results(dat, name=nameTag), stringsAsFactors = FALSE, check.names = FALSE)
  res1$Ensembl <- rownames(res1)
  res1[which(res1$log2FoldChange >= log2(fc) & res1$padj < pthreshold),'sig'] <- 'up'
  res1[which(res1$log2FoldChange <= -log2(fc) & res1$padj < pthreshold),'sig'] <- 'down'
  res1[which(abs(res1$log2FoldChange) <= log2(fc) | res1$padj >= pthreshold),'sig'] <- 'none'
  res1_select <- subset(res1, sig %in% c('up', 'down'))
  res <- merge(res1_select, expData, by = "Ensembl", all.x = T)
  return(res)
}




res_NT_Mock <- filterResult(mm$DEG, "condition_NT_vs_Mock", 2, 0.01)

res_T_Mock <- filterResult(mm$DEG, "condition_T_vs_Mock", 2, 0.01)




# Mock vs NT

