###
 # @Descripttion: Data Analysis
 # @Author: Qun Li
 # @Email: qun.li@ki.se
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

# load tpm, matchedID and sampleinfo
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

# round down expected Count was used to perform DEG analysis
# https://pubmed.ncbi.nlm.nih.gov/35256454/

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
qunDEGAnalysis <- function(dat, condition_CT, num_CT, condition_Treat, num_Treat, fc, pthreshold){
  coldata <- data.frame(condition = factor(c(rep(condition_CT, num_CT), rep(condition_Treat, num_Treat))), levels = c(rep(condition_CT, num_CT), rep(condition_Treat, num_Treat)))
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
  return(res1_select)
}

# expression filter
# exculde genes with fpkm <= 1 



#
expData_test <- round(expData[,21:26])
rownames(expData_test) <- expData$Ensembl
dat = expData_test
condition_CT = "Mock"
num_CT = 3
condition_Treat = "NT"
num_Treat = 3
fc = 2
pthreshold = 0.01

qunDEGAnalysis(expData_test, "Mock", 3, "NT", 3, 2, 0.01)
# PCA Analysis


# DEG for infect and uninfect cells
