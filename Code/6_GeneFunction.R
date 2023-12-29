###
# @Descripttion: function of DE genes
# @Author: LiQun
# @Email: qun.li@ki.se
# @Date: 29 Dec 2023 (10:31:45)
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

# load libraries
{
  library(ggplot2)
  library(clusterProfiler)
  source("./Code/basicFunction.R")
}

# load data
{
  # Ensembl-gene pair
  Ensembl_gene <- read.table("./Data/SourceData/matchedID.txt", header = T, sep = "\t")
  
  # exp this study
  exp_total <- read.table("./Data/ExpRNAseq/exp_Data_merge.txt", header = T, sep = "\t")
  exp_fpkm <- exp_total[,c(1,11:20)]
  
  # exp public data
  # RNAseq analysis of PBMCs in a group of 17 COVID-19 subjects and 17 healthy controls
  exp_GSE152418 <- read.table("./Data/SourceData/HumanPublicCovid19/GSE152418_p20047_Study1_RawCounts.txt", header = T, sep = "\t")
  colnames(exp_GSE152418)[1] <- "Ensembl"
  
  # 126 samples were analyzed in total with 100 COVID-19 patients and 26 non-COVID-19
  exp_GSE157103 <- read.table("./Data/SourceData/HumanPublicCovid19/GSE157103_genes.tpm.tsv", header = T, sep = "\t")
  colnames(exp_GSE157103)[1] <- "Gene"
  
  # Bulk RNA-seq of choroid plexus organoids (CPOs) was performed on mock 72 hours post-infection (hpi), SARS-CoV-2 24 hpi, and SARS-CoV-2 72 hpi samples.
  exp_GSE157852 <- read.table("./Data/SourceData/HumanPublicCovid19/GSE157852_CPO_RawCounts.txt", header = T, sep = "\t")
  
  # Autopsy samples from patients deceased due to SARS-Cov2 infection were collected for total RNA-seq analysis to assess viral load and immune response.
  # lung/heart/liver/kidney/bowel/marrow
  exp_GSE150316 <- read.table("./Data/SourceData/HumanPublicCovid19/GSE150316_RawCounts_Final.txt", header = T, sep = "\t")
  
  # go
  go_mock_NT_up <- read.table("./Results/Table/GO/GO_Data_Mock_NT_Up_strict.txt", header = T, sep = "\t")
  go_mock_NT_down <- read.table("./Results/Table/GO/GO_Data_Mock_NT_Down_strict.txt", header = T, sep = "\t")
  
  go_NT_T_up <- read.table("./Results/Table/GO/GO_Data_NT_T_Up.txt", header = T, sep = "\t")
  go_NT_T_down <- read.table("./Results/Table/GO/GO_Data_NT_T_Down.txt", header = T, sep = "\t")
}

# plot total Go
{
  p_go_mock_NT_up <- qunplotGO(go_mock_NT_up[1:30,])
  p_go_mock_NT_down <- qunplotGO(go_mock_NT_down[1:30,])
  ggsave(filename = "./Results/Figure/6_go_plot_Mock_NT_up_total.pdf", p_go_mock_NT_up, width = 10.29, height = 6.56)
  ggsave(filename = "./Results/Figure/6_go_plot_Mock_NT_down_total.pdf", p_go_mock_NT_down, width = 10.29, height = 6.56)
  
  p_go_mock_NT_up <- qunplotGO(go_mock_NT_up[1:30,])
  p_go_mock_NT_down <- qunplotGO(go_mock_NT_down[1:30,])
  ggsave(filename = "./Results/Figure/6_go_plot_Mock_NT_up_total.pdf", p_go_mock_NT_up, width = 10.29, height = 6.56)
  ggsave(filename = "./Results/Figure/6_go_plot_Mock_NT_down_total.pdf", p_go_mock_NT_down, width = 10.29, height = 6.56)
}

# plot selective go
{
  
}