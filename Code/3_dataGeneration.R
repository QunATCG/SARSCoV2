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
  exp_expectedcount <- read.table("./Data/ExpRNAseq/Human_Covid19_removerRNA.gene.expectedCount", header = T, sep = "\t", stringsAsFactors = F)
  exp_tpm <- read.table("./Data/ExpRNAseq/Human_Covid19_removerRNA.gene.tpm", header = T, sep = "\t", stringsAsFactors = F)
  exp_fpkm <- read.table("./Data/ExpRNAseq/Human_Covid19_removerRNA.gene.fpkm",header = T, sep = "\t", stringsAsFactors = F)
  exp_colnames_common <- c("Mock_1", "Mock_2", "Mock_3", "NT_1", "NT_2", "NT_3", "T_1", "T_2", "T_3")
  colnames(exp_expectedcount) <- c("Ensembl",paste(exp_colnames_common, "expectedreadcount", sep = "_"))
  colnames(exp_tpm) <- c("Ensembl",paste(exp_colnames_common, "tpm", sep = "_"))
  colnames(exp_fpkm) <- c("Ensembl",paste(exp_colnames_common, "fpkm", sep = "_"))
  exp_Ensembl_symbl_readcount <- merge(exp_expectedcount, ensembl_symbol, by = "Ensembl", all.x = T)
  exp_Ensembl_symbl_readcount_tpm <- merge(exp_Ensembl_symbl_readcount, exp_tpm, by = "Ensembl", all.x = T)
  exp_Ensembl_symbl_readcount_tpm_fpkm <- merge(exp_Ensembl_symbl_readcount_tpm, exp_fpkm, by = "Ensembl", all.x = T)
  write.table(exp_Ensembl_symbl_readcount_tpm_fpkm, "./Data/ExpRNAseq/exp_Data_merge.txt", quote = F, sep = "\t", row.names = F)
}
