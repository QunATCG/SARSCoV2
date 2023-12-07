###
# @Descripttion: Data Generation
# @Author: Qun Li
# @Email: qun.li@ki.se/liqun95@163.com
# @Date: 2023-11-23 14:21:57
### 

# set working env
rm(list = ls())
# set your working path 
setwd("/Users/liqun/Desktop/Projects/Covid19/Data/Code/SARSCoV2/")
# check your working path
dirNow <- getwd()

#### generate data from raw files
#### you can skip this generation step
{
  #### STAR_RSEM pipeline
  # load tpm/fpkm/expectedCount/readCount
  # tpm/fpkm/expectedCount were from RSEM (https://github.com/deweylab/RSEM)
  # readCount was from featureCounts (https://subread.sourceforge.net/featureCounts.html)
  # readCount or round down expected Count was used to perform DEG analysis (https://pubmed.ncbi.nlm.nih.gov/35256454/)
  # total
  expTPM_Total_Star <- read.table("./Data/ExpRNAseq/STAR/raw/total.gene.tpm", header = T,
                                  sep = "\t", stringsAsFactors = FALSE)
  expFPKM_Total_Star <- read.table("./Data/ExpRNAseq/STAR/raw/total.gene.fpkm", header = T,
                                   sep = "\t", stringsAsFactors = FALSE)
  expExpextedCount_Total_Star <- read.table("./Data/ExpRNAseq/STAR/raw/total.gene.expected.count", header = T,
                                            sep = "\t", stringsAsFactors = FALSE)
  expReadCount_Total_Star <- read.table("./Data/ExpRNAseq/STAR/raw/total.counts.geneID.out.txt", header = T,
                                        sep = "\t", stringsAsFactors = FALSE)
  
  # covid19
  expTPM_Covid19_Star <- read.table("./Data/ExpRNAseq/STAR/raw/Covid19_Only.gene.tpm", header = T,
                                    sep = "\t", stringsAsFactors = FALSE)
  expFPKM_Covid19_Star <- read.table("./Data/ExpRNAseq/STAR/raw/Covid19_Only.gene.fpkm", header = T,
                                     sep = "\t", stringsAsFactors = FALSE)
  expExpextedCount_Covid19_Star <- read.table("./Data/ExpRNAseq/STAR/raw/Covid19_Only.gene.count", header = T,
                                              sep = "\t", stringsAsFactors = FALSE)
  
  
  #### Hisat2_Stringtie pipeline
  # load tpm/fpkm/expectedCount/readCount
  # tpm/fpkm/readCount were from stringtie (https://ccb.jhu.edu/software/stringtie/)
  # total
  expInTotal <- read.table("./Data/ExpRNAseq/Hisat2/raw/exp_all_sample_covid_total.tsv", header = T,
                           sep = "\t", stringsAsFactors = FALSE)
  expTPM_Total_Hisat2 <- expInTotal[,c(1,grep("TPM",colnames(expInTotal)))]
  expFPKM_Total_Hisat2 <- expInTotal[,c(1,grep("FPKM",colnames(expInTotal)))]
  expReadCount_Total_Hisat2 <- read.table("./Data/ExpRNAseq/Hisat2/raw/total_gene_readcount_new.csv", header = T,
                                          sep = "\t", stringsAsFactors = FALSE)[,c(1,3:11)]
  
  # Covid19
  expInCovid19 <- read.table("./Data/ExpRNAseq/Hisat2/raw/exp_all_sample_covid_only.tsv", header = T,
                             sep = "\t", stringsAsFactors = FALSE)
  expTPM_Covid19_Hisat2 <- expInCovid19[,c(1,grep("TPM",colnames(expInCovid19)))]
  expFPKM_Covid19_Hisat2 <- expInCovid19[,c(1,grep("FPKM",colnames(expInCovid19)))]
  expReadCount_Covid19_Hisat2 <- read.table("./Data/ExpRNAseq/Hisat2/raw/Covid19_gene_readcount_new.csv", header = T,
                                            sep = "\t", stringsAsFactors = FALSE)[,c(1,3:11)]
  
  # matchedID and sampleinfo
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
  
  # rename the colnames of exp
  colNamesTPM <- c("Ensembl",paste(c(sampleInfo$InfoTag),"TPM", sep = "_"))
  colNamesFPKM <- c("Ensembl",paste(c(sampleInfo$InfoTag),"FPKM", sep = "_"))
  colNamesExpectedCount <- c("Ensembl",paste(c(sampleInfo$InfoTag),"ExpectedCount", sep = "_"))
  colNamesReadCount <- c("Ensembl",paste(c(sampleInfo$InfoTag),"ReadCount", sep = "_"))
  
  colnames(expTPM_Total_Star) <- colNamesTPM
  colnames(expTPM_Covid19_Star) <- colNamesTPM
  colnames(expTPM_Total_Hisat2) <- colNamesTPM
  colnames(expTPM_Covid19_Hisat2) <- colNamesTPM
  
  colnames(expFPKM_Total_Star) <- colNamesFPKM
  colnames(expFPKM_Covid19_Star) <- colNamesFPKM
  colnames(expFPKM_Total_Hisat2) <- colNamesFPKM
  colnames(expFPKM_Covid19_Hisat2) <- colNamesFPKM
  
  colnames(expExpextedCount_Total_Star) <- colNamesExpectedCount
  colnames(expExpextedCount_Covid19_Star) <- colNamesExpectedCount
  
  colnames(expReadCount_Total_Star) <- colNamesReadCount
  #colnames(expReadCount_Covid19_Star) <- colNamesReadCount
  colnames(expReadCount_Total_Hisat2) <- colNamesReadCount
  colnames(expReadCount_Covid19_Hisat2) <- colNamesReadCount
  
  
  # match Ensembl and Gene
  expDataStar_Total_TPM <- merge(expTPM_Total_Star, matchedID, by = "Ensembl", all.x = T)
  expDataStar_Total_TPM_FPKM <- merge(expDataStar_Total_TPM, expFPKM_Total_Star, by = "Ensembl", all.x = T)
  expDataStar_Total_TPM_FPKM_ExpectedCount <- merge(expDataStar_Total_TPM_FPKM, expExpextedCount_Total_Star, by = "Ensembl", all.x = T)
  expDataStar_Total_TPM_FPKM_ExpectedCount_ReadCount <- merge(expDataStar_Total_TPM_FPKM_ExpectedCount, expReadCount_Total_Star, by = "Ensembl", all.x = T)
  expDataStar_Total <- expDataStar_Total_TPM_FPKM_ExpectedCount_ReadCount
  
  expDataStar_Covid19_TPM <- merge(expTPM_Covid19_Star, matchedID, by = "Ensembl", all.x = T)
  expDataStar_Covid19_TPM_FPKM <- merge(expDataStar_Covid19_TPM, expFPKM_Covid19_Star, by = "Ensembl", all.x = T)
  expDataStar_Covid19_TPM_FPKM_ExpectedCount <- merge(expDataStar_Covid19_TPM_FPKM, expExpextedCount_Covid19_Star, by = "Ensembl", all.x = T)
  expDataStar_Covid19 <- expDataStar_Covid19_TPM_FPKM_ExpectedCount
  
  
  expDataHisat2_Total_TPM <- merge(expTPM_Total_Hisat2, matchedID, by = "Ensembl", all.x = T)
  expDataHisat2_Total_TPM_FPKM <- merge(expDataHisat2_Total_TPM, expFPKM_Total_Hisat2, by = "Ensembl", all.x = T)
  expDataHisat2_Total_TPM_FPKM_ReadCount <- merge(expDataHisat2_Total_TPM_FPKM, expReadCount_Total_Hisat2, by = "Ensembl", all.x = T)
  expDataHisat2_Total <- expDataHisat2_Total_TPM_FPKM_ReadCount
  
  expDataHisat2_Covid19_TPM <- merge(expTPM_Covid19_Hisat2, matchedID, by = "Ensembl", all.x = T)
  expDataHisat2_Covid19_TPM_FPKM <- merge(expDataHisat2_Covid19_TPM, expFPKM_Covid19_Hisat2, by = "Ensembl", all.x = T)
  expDataHisat2_Covid19_TPM_FPKM_ReadCount <- merge(expDataHisat2_Covid19_TPM_FPKM, expReadCount_Covid19_Hisat2, by = "Ensembl", all.x = T)
  expDataHisat2_Covid19 <- expDataHisat2_Covid19_TPM_FPKM_ReadCount
  
  # save data
  write.table(expDataStar_Total, "./Data/ExpRNAseq/STAR/Star_exp_total.txt", sep = "\t", row.names = F, quote = FALSE)
  write.table(expDataStar_Covid19, "./Data/ExpRNAseq/STAR/Star_exp_Covid19.txt", sep = "\t", row.names = F, quote = FALSE)
  write.table(expDataHisat2_Total, "./Data/ExpRNAseq/Hisat2/Hisat2_exp_total.txt", sep = "\t", row.names = F, quote = FALSE)
  write.table(expDataHisat2_Covid19, "./Data/ExpRNAseq/Hisat2/Hisat2_exp_Covid19.txt", sep = "\t", row.names = F, quote = FALSE)
}
