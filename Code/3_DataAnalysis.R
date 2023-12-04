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
expTPM <- read.table("./Data/ExpRNAseq/Covid19.gene.tpm", header = T,
                    sep = "\t", stringsAsFactors = FALSE)
expCount <- read.table("./Data/ExpRNAseq/Covid19.gene.count", header = T,
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
colNamesExp <- c("Ensembl",paste(c(sampleInfo$InfoTag),"TPM", sep = "_"))
colNamesCount <- c("Ensembl",paste(c(sampleInfo$InfoTag),"Count", sep = "_"))
colnames(expTPM) <- colNamesExp
colnames(expCount) <- colNamesCount

# match Ensembl and Gene
expDataTPM <- merge(expTPM, matchedID, by = "Ensembl", all.x = T)
expDataTPMCount <- merge(expDataTPM, expCount, by = "Ensembl", all.x = T)
expData <- expDataTPMCount

# PCA Analysis


# DEG for infect and uninfect cells
