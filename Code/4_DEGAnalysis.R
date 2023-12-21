###
 # @Descripttion: Data Analysis
 # @Author: Qun Li
 # @Email: qun.li@ki.se/liqun95@163.com
 # @Date: 2023-11-23 14:21:57
### 

# DEG analysis was performed using R package DESeq2
# ploting was performed using R package pheatmap

# env setting
{
  # set working env
  rm(list = ls())
  # set your working path 
  setwd("/Users/liqun/Desktop/Projects/Covid19/Data/Code/SARSCoV2/")
  # check your working path
  dirNow <- getwd()
}

# load libraries
{
  library(DESeq2)
  library(pheatmap)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  source("./Code/basicFunction.R")
}

# load exp data
{
  exp_data_total <- read.table("./Data/ExpRNAseq/exp_Data_merge.txt", header = T, sep = "\t", stringsAsFactors = F, row.names = 1)
  exp_data_readcounts <- exp_data_total[,1:9]
  # remove all 0 read counts genes
  exp_data_readcounts <- exp_data_readcounts[rowSums(exp_data_readcounts) !=0, ]
}

# PCA analysis
# get normalized read counts for PCA plotting
{
  data_PCA <- t(exp_data_readcounts)
  # remove 0
  data_PCA_prcomp <- prcomp(data_PCA, scale=TRUE)
  # calculate percentage of PCs
  data_PCA_var <- data_PCA_prcomp$sdev^2
  data_PCA_per <- round(data_PCA_var/sum(data_PCA_var)*100,1)
  write.table(data_PCA_per, "./Data/ExpRNAseq/exp_Data_PCA_percentage.txt", quote = F, sep = "\t", row.names = T)
  # PCs 
  data_PCA_output <- data.frame(data_PCA_prcomp$x)
  write.table(data_PCA_output, "./Data/ExpRNAseq/exp_Data_PCA.txt", quote = F, sep = "\t", row.names = T)
  
  # correlation of samples
  data_sample_cor <- cor(exp_data_readcounts)
  write.table(data_sample_cor, "./Data/ExpRNAseq/exp_Data_Correlation.txt", quote = F, sep = "\t", row.names = T)
}

# DE analysis
{
  # DE
  DE_Mock_NT <- qunDEGAnalysisTwoCondition(dat = exp_data_readcounts[,1:6], condition_CT = "Mock", num_CT = 3, condition_Treat = "NT", num_Treat = 3)
  DE_Mock_T  <- qunDEGAnalysisTwoCondition(dat = exp_data_readcounts[,c(1:3,7:9)], condition_CT = "Mock", num_CT = 3, condition_Treat = "T", num_Treat = 3)
  DE_NT_T    <- qunDEGAnalysisTwoCondition(dat = exp_data_readcounts[,4:9], condition_CT = "NT", num_CT = 3, condition_Treat = "T", num_Treat = 3)
  
  # mark results
  DE_Mock_NT_result <- markResult(dat = DE_Mock_NT$DEG, fc = 1.5, pthreshold = 0.05, tagItem = "treatVScontrol")
  DE_Mock_T_result  <- markResult(dat = DE_Mock_T$DEG,  fc = 1.5, pthreshold = 0.05, tagItem = "treatVScontrol")
  DE_NT_T_result    <- markResult(dat = DE_NT_T$DEG,    fc = 1.5, pthreshold = 0.05, tagItem = "treatVScontrol")
  
  DE_Mock_NT_result$Ensembl <- rownames(DE_Mock_NT_result)
  DE_Mock_T_result$Ensembl <- rownames(DE_Mock_T_result)
  DE_NT_T_result$Ensembl <- rownames(DE_NT_T_result)
  
  exp_data_total$Ensembl <- rownames(exp_data_total)
  
  # merge exp files
  DE_Mock_NT_result_exp <- merge(DE_Mock_NT_result, exp_data_total, by = "Ensembl", all.x = T)
  DE_Mock_T_result_exp <-  merge(DE_Mock_T_result,  exp_data_total, by = "Ensembl", all.x = T)
  DE_NT_T_result_exp <-    merge(DE_NT_T_result,    exp_data_total, by = "Ensembl", all.x = T)
  
  # exclude baseMean = 0
  DE_Mock_NT_result_exp <- DE_Mock_NT_result_exp[DE_Mock_NT_result_exp$baseMean != 0,]
  DE_Mock_T_result_exp  <- DE_Mock_T_result_exp[DE_Mock_T_result_exp$baseMean != 0,]
  DE_NT_T_result_exp    <- DE_NT_T_result_exp[DE_NT_T_result_exp$baseMean != 0,]
  
  # Output
  write.table(DE_Mock_NT_result_exp, "./Data/ExpRNAseq/DE_Data_Mock_NT.txt", quote = F, sep = "\t", row.names = F)
  write.table(DE_Mock_T_result_exp, "./Data/ExpRNAseq/DE_Data_Mock_T.txt", quote = F, sep = "\t", row.names = F)
  write.table(DE_NT_T_result_exp, "./Data/ExpRNAseq/DE_Data_NT_T.txt", quote = F, sep = "\t", row.names = F)
}

# GO analysis
{
  # GO
  GO_Mock_NT_Up   <- qunGO(dat = DE_Mock_NT_result_exp, "up")
  GO_Mock_NT_Down <- qunGO(dat = DE_Mock_NT_result_exp, "down")
  
  GO_Mock_T_Up    <- qunGO(dat = DE_Mock_T_result_exp, "up")
  GO_Mock_T_Down  <- qunGO(dat = DE_Mock_T_result_exp, "down")
  
  GO_NT_T_Up      <- qunGO(dat = DE_NT_T_result_exp, "up")
  GO_NT_T_Down    <- qunGO(dat = DE_NT_T_result_exp, "down")
  
  # Output
  write.table(GO_Mock_NT_Up,   "./Data/ExpRNAseq/GO_Data_Mock_NT_Up.txt", quote = F, sep = "\t", row.names = F)
  write.table(GO_Mock_NT_Down, "./Data/ExpRNAseq/GO_Data_Mock_NT_Down.txt", quote = F, sep = "\t", row.names = F)
  
  write.table(GO_Mock_T_Up,    "./Data/ExpRNAseq/GO_Data_Mock_T_Up.txt", quote = F, sep = "\t", row.names = F)
  write.table(GO_Mock_T_Down,  "./Data/ExpRNAseq/GO_Data_Mock_T_Down.txt", quote = F, sep = "\t", row.names = F)
  
  write.table(GO_NT_T_Up,      "./Data/ExpRNAseq/GO_Data_NT_T_Up.txt", quote = F, sep = "\t", row.names = F)
  write.table(GO_NT_T_Down,    "./Data/ExpRNAseq/GO_Data_NT_T_Down.txt", quote = F, sep = "\t", row.names = F)
}