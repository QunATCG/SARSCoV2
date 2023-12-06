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

# load tpm/fpkm/expectedCount/readCount, matchedID and sampleinfo
# tpm/fpkm/expectedCount were from RSEM (https://github.com/deweylab/RSEM)
# readCount was from featureCounts (https://subread.sourceforge.net/featureCounts.html)
# readCount or round down expected Count was used to perform DEG analysis (https://pubmed.ncbi.nlm.nih.gov/35256454/)

# 
expTPM <- read.table("./Data/ExpRNAseq/STAR/Covid19.gene.tpm", header = T,
                    sep = "\t", stringsAsFactors = FALSE)
expFPKM <- read.table("./Data/ExpRNAseq/STAR/Covid19.gene.fpkm", header = T,
                      sep = "\t", stringsAsFactors = FALSE)
expExpextedCount <- read.table("./Data/ExpRNAseq/STAR/Covid19.gene.count", header = T,
                       sep = "\t", stringsAsFactors = FALSE)
expReadCount <- read.table("./Data/ExpRNAseq/STAR/Covid19.counts.geneID.out.txt", header = T,
                           sep = "\t", stringsAsFactors = FALSE)
matchedID <- read.table("./Data/matchedID.txt", header = T, sep = "\t",
                        stringsAsFactors = FALSE)
sampleInfo <- read.table("./Data/Sample_information.txt", header = T,
                        sep = "\t", stringsAsFactors = FALSE)



# process NA items of matchedID
for (i in 1:nrow(matchedID)){
  if (matchedID[i,]$Gene == ""){
    matchedID[i,]$Gene = matchedID[i,]$Ensembl
  }
}

# rename the colnames of expTPM
colNamesTPM <- c("Ensembl",paste(c(sampleInfo$InfoTag),"TPM", sep = "_"))
colNamesFPKM <- c("Ensembl",paste(c(sampleInfo$InfoTag),"FPKM", sep = "_"))
colNamesExpectedCount <- c("Ensembl",paste(c(sampleInfo$InfoTag),"ExpectedCount", sep = "_"))
colNamesReadCount <- c("Ensembl",paste(c(sampleInfo$InfoTag),"ReadCount", sep = "_"))
colnames(expTPM) <- colNamesTPM
colnames(expFPKM) <- colNamesFPKM
colnames(expExpextedCount) <- colNamesExpectedCount
colnames(expReadCount) <- colNamesReadCount

# match Ensembl and Gene
expDataTPM <- merge(expTPM, matchedID, by = "Ensembl", all.x = T)
expDataTPM_FPKM <- merge(expDataTPM, expFPKM, by = "Ensembl", all.x = T)
expDataTPM_FPKM_ExpectedCount <- merge(expDataTPM_FPKM, expExpextedCount, by = "Ensembl", all.x = T)
expDataTPM_FPKM_ExpectedCount_ReadCount <- merge(expDataTPM_FPKM_ExpectedCount, expReadCount, by = "Ensembl", all.x = T)
expData <- expDataTPM_FPKM_ExpectedCount_ReadCount

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

# expression filter
# exculde genes with fpkm <= 1 
expDataOver1 <- expData[which(rowMeans(expData[,c("Mock_1_FPKM", "Mock_2_FPKM", "Mock_3_FPKM")]) >=1 | 
                      rowMeans(expData[,c("T_1_FPKM", "T_2_FPKM", "T_3_FPKM")]) >=1 | 
                      rowMeans(expData[,c("NT_1_FPKM", "NT_2_FPKM", "NT_3_FPKM")]) >=1),]

# PCA
expData_PCA <- expDataOver1[,30:38]
rownames(expData_PCA) <- expDataOver1$Ensembl
mm <- qunDEGAnalysisThreeCondition(round(expData_PCA), "Mock", 3, "NT", 3, "T", 3, 2, 0.01)
plotDispEsts(mm$DEG, main="Dispersion plot")
# get normalized readcounts from DESeq2 results
rld <- rlog(mm$DEG)
plotPCA(rld, ntop = 3000)

rld_mat <- assay(rld)
rld_cor <- cor(rld_mat) 
pheatmap(rld_cor)


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

