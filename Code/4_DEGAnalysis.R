###
 # @Descripttion: Data Analysis
 # @Author: Qun Li
 # @Email: qun.li@ki.se/liqun95@163.com
 # @Date: 2023-11-23 14:21:57
### 

{
  # DEG analysis was performed using R package DESeq2
  # ploting was performed using R package pheatmap
  # <RSEM-DESeq2> this method was used by some groups, such as
  # https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
  # https://pubmed.ncbi.nlm.nih.gov/38097578/
  # https://pubmed.ncbi.nlm.nih.gov/38092806/
  # https://pubmed.ncbi.nlm.nih.gov/33870146/
  # https://pubmed.ncbi.nlm.nih.gov/34528097/
  # https://pubmed.ncbi.nlm.nih.gov/38049398/
}

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
  library(tximport)
  library(pheatmap)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  source("./Code/basicFunction.R")
}

# load exp data
{
  Ensembl_gene <- read.table("./Data/SourceData/matchedID.txt", header = T, sep = "\t")
  exp_data_total <- read.table("./Data/ExpRNAseq/exp_Data_merge.txt", header = T, sep = "\t", stringsAsFactors = F, row.names = 1)
  exp_data_readcounts <- round(exp_data_total[,1:9])
  exp_data_fpkm <- exp_data_total[,11:19]
  exp_data_readcounts_total <- read.table("./Data/ExpRNAseq/finalExp/Human_Covid19_removerRNA.gene.readCounts", header = T, sep = "\t", row.names = 1)
  # calculate variance
  exp_data_readcounts$Variance <- apply(exp_data_readcounts, 1, var)
  exp_data_fpkm$Variance <- apply(exp_data_fpkm, 1, var)
  # check na rows
  exp_data_readcounts[is.na(exp_data_readcounts),]
  exp_data_fpkm[is.na(exp_data_fpkm),]
  # sars genes
  sars_genes <- Ensembl_gene$Ensembl[grep("ENSSASG", Ensembl_gene$Ensembl)]
}

# PCA analysis
# get normalized read counts for PCA plotting
{
  ### based on read counts
  {
    exp_data_readcounts_PCA <- exp_data_readcounts_total
    colnames(exp_data_readcounts_PCA) <- c("Mock_1", "Mock_2", "Mock_3", "NT_1", "NT_2", "NT_3", "T_1","T_2", "T_3")
    exp_data_readcounts_PCA$Variance <- apply(exp_data_readcounts_PCA, 1, var)
    # exclude sars genes
    # exp_data_readcounts_PCA <- exp_data_readcounts[!rownames(exp_data_readcounts) %in% sars_genes,]
    # remove constant/zero reads genes
    data_readcount_PCA <- as.data.frame(t(exp_data_readcounts_PCA[exp_data_readcounts_PCA$Variance != 0,][,1:9]))
    data_readcount_PCA_prcomp <- prcomp(data_readcount_PCA, scale. = T)

    # calculate percentage of PCs
    data_readcount_PCA_var <- data_readcount_PCA_prcomp$sdev^2
    data_readcount_PCA_per <- round(data_readcount_PCA_var/sum(data_readcount_PCA_var)*100,1)
    write.table(data_readcount_PCA_per, "./Results/Table/exp_Data_readcount_PCA_percentage.txt", quote = F, sep = "\t", row.names = T)
    # PCs
    data_readcount_PCA_output <- data.frame(data_readcount_PCA_prcomp$x)
    # Check PCA plot
    plot(data_readcount_PCA_output[,1], data_readcount_PCA_output[,2])
    write.table(data_readcount_PCA_output, "./Results/Table/exp_Data_readcount_PCA.txt", quote = F, sep = "\t", row.names = T)

    # correlation
    exp_data_readcounts_cor <- cor(exp_data_readcounts_PCA[,1:9])
    pheatmap(exp_data_readcounts_cor)
    write.table(exp_data_readcounts_cor, "./Results/Table/exp_Data_readcount_cor.txt", quote = F, sep = "\t", row.names = T)
  }
  
  # ### based on fpkm
  # {
  #   # exclude sars genes
  #   exp_data_fpkm_PCA <- exp_data_fpkm[!rownames(exp_data_fpkm) %in% sars_genes,]
  #   exp_data_fpkm_PCA <- na.omit(exp_data_fpkm_PCA)
  #   # remove constant/zero reads genes
  #   data_fpkm_PCA <- as.data.frame(t(log2(exp_data_fpkm_PCA[exp_data_fpkm_PCA$Variance != 0,][,1:9] + 1)))
  #   data_fpkm_PCA <- data_fpkm_PCA[,apply(data_fpkm_PCA,2,var) != 0]
  #   data_fpkm_PCA_prcomp <- prcomp(data_fpkm_PCA, scale. = T)
  # 
  #   # calculate percentage of PCs
  #   data_fpkm_PCA_var <- data_fpkm_PCA_prcomp$sdev^2
  #   data_fpkm_PCA_per <- round(data_fpkm_PCA_var/sum(data_fpkm_PCA_var)*100,1)
  #   write.table(data_fpkm_PCA_per, "./Results/Table/exp_Data_fpkm_PCA_percentage.txt", quote = F, sep = "\t", row.names = T)
  #   # PCs
  #   data_fpkm_PCA_output <- data.frame(data_fpkm_PCA_prcomp$x)
  #   # Check PCA plot
  #   plot(data_fpkm_PCA_output[,1], data_fpkm_PCA_output[,2])
  #   write.table(data_fpkm_PCA_output, "./Results/Table/exp_Data_fpkm_PCA.txt", quote = F, sep = "\t", row.names = T)
  # 
  #   # correlation
  #   exp_data_fpkm_cor <- cor(exp_data_fpkm[,c(1:9)])
  #   pheatmap(exp_data_fpkm_cor, clustering_method = "average")
  #   write.table(exp_data_fpkm_cor, "./Results/Table/exp_Data_fpkm_cor.txt", quote = F, sep = "\t", row.names = T)
  # }
}

# DE analysis
{
  # DE
  # featureCount-DESeq2
  {
    # sample info
    samples_info <- read.table("./Data/SourceData/Sample_information.txt", header = T, sep = "\t")
    samples_info$treatment <- c(rep("Mock", 3), rep("NT", 3), rep("T",3))
    rownames(samples_info) <- samples_info$Rowname
    
    # DE of Mock and NT (exclude covid read counts)
    # load data
    exp_data_readcounts_DE_Mock_NT <- exp_data_readcounts[exp_data_readcounts$Variance != 0,][,1:6]
    
    # DE
    data_dds_Mock_NT <- DESeqDataSetFromMatrix(countData = exp_data_readcounts_DE_Mock_NT, colData = samples_info[1:6,], design= ~treatment)
    data_dds_Mock_NT <- data_dds_Mock_NT[rowSums(counts(data_dds_Mock_NT)) > 1,]
    data_dds_DESeq_Mock_NT <- DESeq(data_dds_Mock_NT, fitType = "mean")
    DE_Mock_NT <- results(data_dds_DESeq_Mock_NT, contrast=c("treatment","Mock","NT"))
    
    
    # DE of NT and T (include covid read counts)
    # load data
    exp_data_readcounts_DE_NT_T <- exp_data_readcounts_PCA[exp_data_readcounts_PCA$Variance != 0,][,4:9]
    colnames(exp_data_readcounts_DE_NT_T) <- c("NT_1_rawcount", "NT_2_rawcount", "NT_3_rawcount",
                                               "T_1_rawcount", "T_2_rawcount", "T_3_rawcount")
    
    # DE
    data_dds_NT_T <- DESeqDataSetFromMatrix(countData = exp_data_readcounts_DE_NT_T, colData = samples_info[4:9,], design= ~treatment)
    
    data_dds_NT_T <- data_dds_NT_T[rowSums(counts(data_dds_NT_T)) > 1,]
    data_dds_DESeq_NT_T <- DESeq(data_dds_NT_T, fitType = "mean")
    DE_NT_T <- results(data_dds_DESeq_NT_T, contrast=c("treatment","NT","T"))
    
    # mark results
    DE_Mock_NT_result <- markResult(dat = data.frame(DE_Mock_NT), fc = 1.5, pthreshold = 0.05, tagItem = "controlVStreat")
    DE_NT_T_result    <- markResult(dat = data.frame(DE_NT_T),    fc = 1.5, pthreshold = 0.05, tagItem = "controlVStreat")
    
    DE_Mock_NT_result$Ensembl <- rownames(DE_Mock_NT_result)
    DE_NT_T_result$Ensembl <- rownames(DE_NT_T_result)
    
    exp_data_total$Ensembl <- rownames(exp_data_total)
    
    # remove NA
    #DE_Mock_NT_result <- DE_Mock_NT_result[!is.na(DE_Mock_NT_result$padj),]
    #DE_Mock_T_result <- DE_Mock_T_result[!is.na(DE_Mock_T_result$padj),]
    #DE_NT_T_result <- DE_NT_T_result[!is.na(DE_NT_T_result$padj),]

    # merge exp files
    DE_Mock_NT_result_exp <- merge(DE_Mock_NT_result, exp_data_total, by = "Ensembl", all.x = T)
    DE_NT_T_result_exp <-    merge(DE_NT_T_result,    exp_data_total, by = "Ensembl", all.x = T)

    # exclude baseMean = 0
    DE_Mock_NT_result_exp <- DE_Mock_NT_result_exp[DE_Mock_NT_result_exp$baseMean != 0,]
    DE_NT_T_result_exp    <- DE_NT_T_result_exp[DE_NT_T_result_exp$baseMean != 0,]
    
    # exclude readcounts <= 5
    DE_Mock_NT_result_exp <- keepOverReadCount(DE_Mock_NT_result_exp, 5)
    DE_NT_T_result_exp    <- keepOverReadCount(DE_NT_T_result_exp, 5)
    
    # exclude fpkm <= 1
    DE_Mock_NT_result_exp <- keepOverExp(DE_Mock_NT_result_exp, 1)
    DE_NT_T_result_exp <- keepOverExp(DE_NT_T_result_exp, 1)
    
    # Output
    write.table(DE_Mock_NT_result_exp, "./Results/Table/DEG/DE_Data_Mock_NT.txt", quote = F, sep = "\t", row.names = F)
    write.table(DE_NT_T_result_exp, "./Results/Table/DEG/DE_Data_NT_T.txt", quote = F, sep = "\t", row.names = F)
  }
  
  # DE
  ## RSEM-DESeq2
  # https://pubmed.ncbi.nlm.nih.gov/38097578/
  # https://pubmed.ncbi.nlm.nih.gov/38092806/
  # https://pubmed.ncbi.nlm.nih.gov/33870146/
  # https://pubmed.ncbi.nlm.nih.gov/34528097/
  # https://pubmed.ncbi.nlm.nih.gov/38049398/
  # {
  #   # load data
  #   samples_info <- read.table("./Data/SourceData/Sample_information.txt", header = T, sep = "\t")
  #   samples_info$treatment <- c(rep("Mock", 3), rep("NT", 3), rep("T",3))
  #   files_path <- paste("./Data/ExpRNAseq/rsem/", samples_info$NumSeq, ".genes.results", sep = "")
  #   txi.rsem <- tximport(files_path, type = "rsem", txIn = FALSE, txOut = FALSE)
  #   txi.rsem$length[txi.rsem$length == 0] <- 1
  #   
  #   # import data to DESeq2
  #   dds <- DESeqDataSetFromTximport(txi.rsem, colData = samples_info, design = ~ treatment)
  #   #dds <- dds[rowSums(counts(dds)) > 1,]
  #   dds_DESeq <- DESeq(dds)
  #   
  #   DE_Mock_NT <- results(dds_DESeq, contrast=c("treatment","Mock","NT"))
  #   DE_Mock_T <- results(dds_DESeq, contrast=c("treatment","Mock","T"))
  #   DE_NT_T <- results(dds_DESeq, contrast=c("treatment","NT","T")) 
  #   
  #   
  #   # mark results
  #   DE_Mock_NT_result <- markResult(dat = data.frame(DE_Mock_NT), fc = 1.5, pthreshold = 0.05, tagItem = "controlVStreat")
  #   DE_Mock_T_result  <- markResult(dat = data.frame(DE_Mock_T),  fc = 1.5, pthreshold = 0.05, tagItem = "controlVStreat")
  #   DE_NT_T_result    <- markResult(dat = data.frame(DE_NT_T),    fc = 1.5, pthreshold = 0.05, tagItem = "controlVStreat")
  #   
  #   DE_Mock_NT_result$Ensembl <- rownames(DE_Mock_NT_result)
  #   DE_Mock_T_result$Ensembl <- rownames(DE_Mock_T_result)
  #   DE_NT_T_result$Ensembl <- rownames(DE_NT_T_result)
  #   
  #   exp_data_total$Ensembl <- rownames(exp_data_total)
  #   
  #   # remove NA
  #   DE_Mock_NT_result <- DE_Mock_NT_result[!is.na(DE_Mock_NT_result$padj),]
  #   DE_Mock_T_result <- DE_Mock_T_result[!is.na(DE_Mock_T_result$padj),]
  #   DE_NT_T_result <- DE_NT_T_result[!is.na(DE_NT_T_result$padj),]
  #   
  #   
  #   # merge exp files
  #   DE_Mock_NT_result_exp <- merge(DE_Mock_NT_result, exp_data_total, by = "Ensembl", all.x = T)
  #   DE_Mock_T_result_exp <-  merge(DE_Mock_T_result,  exp_data_total, by = "Ensembl", all.x = T)
  #   DE_NT_T_result_exp <-    merge(DE_NT_T_result,    exp_data_total, by = "Ensembl", all.x = T)
  #   
  #   # exclude baseMean = 0
  #   DE_Mock_NT_result_exp <- DE_Mock_NT_result_exp[DE_Mock_NT_result_exp$baseMean != 0,]
  #   DE_Mock_T_result_exp  <- DE_Mock_T_result_exp[DE_Mock_T_result_exp$baseMean != 0,]
  #   DE_NT_T_result_exp    <- DE_NT_T_result_exp[DE_NT_T_result_exp$baseMean != 0,]
  #   
  #   # Output
  #   write.table(DE_Mock_NT_result_exp, "./Results/Table/DEG/DE_Data_Mock_NT.txt", quote = F, sep = "\t", row.names = F)
  #   write.table(DE_Mock_T_result_exp, "./Results/Table/DEG/DE_Data_Mock_T.txt", quote = F, sep = "\t", row.names = F)
  #   write.table(DE_NT_T_result_exp, "./Results/Table/DEG/DE_Data_NT_T.txt", quote = F, sep = "\t", row.names = F)
  # }
}

# GO analysis
{
  # GO
  # For Mock_NT_up genes, extract pvalue <= 0.01 and abs(log2FC) >= 1 genes
  DE_Mock_NT_result_exp_strict <- DE_Mock_NT_result_exp[abs(DE_Mock_NT_result_exp$log2FoldChange) >= 1 & DE_Mock_NT_result_exp$pvalue <= 0.01,]
  GO_Mock_NT_Up   <- qunGO(dat = DE_Mock_NT_result_exp_strict, "up")
  GO_Mock_NT_Down <- qunGO(dat = DE_Mock_NT_result_exp_strict, "down")
  
  GO_NT_T_Up      <- qunGO(dat = DE_NT_T_result_exp, "up")
  GO_NT_T_Down    <- qunGO(dat = DE_NT_T_result_exp, "down")
  
  # Output
  write.table(GO_Mock_NT_Up,   "./Results/Table/GO/GO_Data_Mock_NT_Up.txt", quote = F, sep = "\t", row.names = F)
  write.table(GO_Mock_NT_Down, "./Results/Table/GO/GO_Data_Mock_NT_Down.txt", quote = F, sep = "\t", row.names = F)
  
  write.table(GO_NT_T_Up,      "./Results/Table/GO/GO_Data_NT_T_Up.txt", quote = F, sep = "\t", row.names = F)
  write.table(GO_NT_T_Down,    "./Results/Table/GO/GO_Data_NT_T_Down.txt", quote = F, sep = "\t", row.names = F)
}
