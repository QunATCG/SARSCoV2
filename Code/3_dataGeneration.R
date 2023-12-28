###
# @Descripttion: Data Generation
# @Author: Qun Li
# @Email: qun.li@ki.se/liqun95@163.com
# @Date: 2023-11-23 14:21:57
### 

# env setting
{
  # set working env
  rm(list = ls())
  # set your working path 
  setwd("/Users/liqun/Desktop/Projects/Covid19/Data/Code/SARSCoV2/")
  # check your working path
  dirNow <- getwd()
}


#### generate data from raw files
{
  ensembl_symbol <- read.table("./Data/SourceData/matchedID.txt", header = T, sep = "\t", stringsAsFactors = F)
  sample_info <- read.table("./Data/SourceData/Sample_information.txt", header = T, sep = "\t", stringsAsFactors = F)
  exp_fpkm <- read.table("./Data/ExpRNAseq/finalExp/gene_expression_fpkm.txt", header = T, sep = "\t", stringsAsFactors = F)
  exp_readcounts <- read.table("./Data/ExpRNAseq/finalExp/Human_Covid19_removerRNA.gene.readCounts", header = T, sep = "\t", stringsAsFactors = F)
  exp_colnames_common <- c("Mock_1", "Mock_2", "Mock_3", "NT_1", "NT_2", "NT_3", "T_1", "T_2", "T_3")
  colnames(exp_fpkm) <- c("Ensembl",paste(exp_colnames_common, "fpkm", sep = "_"))
  colnames(exp_readcounts) <- c("Ensembl", paste(exp_colnames_common, "rawcount", sep = "_"))
  exp_Ensembl_symbl_readcount <- merge(exp_readcounts, ensembl_symbol, by = "Ensembl", all.x = T)
  # only keep both
  exp_Ensembl_symbl_readcount_fpkm <- merge(exp_Ensembl_symbl_readcount, exp_fpkm, by = "Ensembl")
  write.table(exp_Ensembl_symbl_readcount_fpkm, "./Data/ExpRNAseq/exp_Data_merge.txt", quote = F, sep = "\t", row.names = F)
}
