###
 # @Descripttion: plotData
 # @Author: Qun Li
 # @Email: qun.li@ki.se/liqun95@163.com
 # @Date: 2023-12-07 14:21:57
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


### load libraries
{
  library(ggplot2)
  library(gplots)
  library(ggrepel)
  library(ggsci)
  library(ggpubr)
  library(pheatmap)
  library(VennDiagram)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  source("./Code/basicFunction.R")
  myColors <- c("grey", "#00A087B2", "#DC0000B2")
}

# load exp data
{
  Ensembl_gene <- read.table("./Data/SourceData/matchedID.txt", header = T, sep = "\t")
  exp_data_total <- read.table("./Data/ExpRNAseq/exp_Data_merge.txt", header = T, sep = "\t", stringsAsFactors = F, row.names = 1)
  fpkm_cols <- c("Mock_1_fpkm", "Mock_2_fpkm", "Mock_3_fpkm", "NT_1_fpkm", 
                 "NT_2_fpkm", "NT_3_fpkm", "T_1_fpkm", "T_2_fpkm", "T_3_fpkm")
  rawcount_cols <- c("Mock_1_rawcount", "Mock_2_rawcount", "Mock_3_rawcount", "NT_1_rawcount", 
                     "NT_2_rawcount", "NT_3_rawcount", "T_1_rawcount", "T_2_rawcount", "T_3_rawcount")
  fpkm_mean_cols <- c("Mock_mean", "NT_mean", "T_mean")
  exp_data_rawcounts <- exp_data_total[,rawcount_cols]
  exp_data_fpkm <- exp_data_total[,fpkm_cols]
  exp_data_fpkm$Mock_mean <- apply(exp_data_fpkm[,1:3], 1, mean)
  exp_data_fpkm$NT_mean <- apply(exp_data_fpkm[,4:6], 1, mean)
  exp_data_fpkm$T_mean <- apply(exp_data_fpkm[,7:9], 1, mean)
  # remove all 0 read counts genes
  exp_data_rawcounts <- exp_data_rawcounts[rowSums(exp_data_rawcounts) !=0, ]
  # SARS genes
  sars_genes <- Ensembl_gene$Ensembl[grep("ENSSASG", Ensembl_gene$Ensembl)]
}

# PCA 
{
  data_PCA <- read.table("./Results/Table/exp_Data_readcount_PCA.txt", header = T, row.names = 1, sep = "\t")
  data_PCA$group <- c(rep("Mock",3), rep("NT",3), rep("T",3))
  data_PCA$label <- c("Mock_1", "Mock_2", "Mock_3", "NT_1", "NT_2", "NT_3", "T_1", "T_2", "T_3")
  data_PCA_percentage <- read.table("./Results/Table/exp_Data_readcount_PCA_percentage.txt", header = T)
  
  p_pca <- ggplot(data = data_PCA, mapping = aes(x = PC1, y = PC2, colour = group)) + geom_point(size = 5) +
    geom_text_repel(aes(label = label), size = 4) +
    scale_color_manual(values = c( Mock = "black", NT = "gray88", T = "darkgrey" )) +
    xlab(paste("PC1: ",data_PCA_percentage$x[1],"% variance", sep = "")) +
    ylab(paste("PC2: ",data_PCA_percentage$x[2],"% variance", sep = "")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  p_pca
  #pdf("./Results/Figure/0_PCA.pdf", width = 7.0, height = 5.36)
  ggsave(filename = "./Results/Figure/0_PCA.pdf", p_pca, width = 7.0, height = 5.36)
  #dev.off()
}

# Correlation
{
  data_cor <- read.table("./Results/Table/exp_Data_readcount_cor.txt", header = T, row.names = 1, sep = "\t")
  p_pheatmap <- pheatmap(data_cor, clustering_method = "median")
  p_pheatmap
  #pdf("./Results/Figure/1_Correlation.pdf", width = 7.0, height = 5.36)
  ggsave(filename = "./Results/Figure/1_Correlation.pdf", p_pheatmap, width = 7.0, height = 5.36)
  #dev.off()
}

# percentage of SARS readcounts
{
  data_rawcounts_total_sars <- read.table("./Data/ExpRNAseq/finalExp/Human_Covid19_removerRNA.gene.readCounts", header = T, sep = "\t", row.names = 1)
  rawcounts_total <- colSums(data_rawcounts_total_sars)
  rawcounts_sars  <- colSums(data_rawcounts_total_sars[rownames(data_rawcounts_total_sars) %in% sars_genes,])
  sars_rawcounts_percentage <- rawcounts_sars/rawcounts_total
  sars_rawcounts_percentage_data <- data.frame(sample = rep(c("Mock", "NT", "T"), each = 3), value = sars_rawcounts_percentage)
  sars_rawcounts_percentage_NT_mean <- mean(sars_rawcounts_percentage_data$value[4:6])
  sars_rawcounts_percentage_T_mean <- mean(sars_rawcounts_percentage_data$value[7:9])
  pvalue_NT_T <- t.test(sars_rawcounts_percentage_data$value[4:6], sars_rawcounts_percentage_data$value[7:9])$p.value
  myComparision1 <- list(c("NT","T"))
  
  p_covidReadCountsPercent <- ggplot(sars_rawcounts_percentage_data,aes(x = sample, y = value, fill = sample)) +
    geom_bar(stat = "summary", fun = mean, width = 0.5) +
    scale_fill_manual(values = c("black", "grey88", "darkgrey")) +
    #stat_summary(mapping = aes(fill = sample),fun = mean, geom = "bar",fun.args = list(mult = 1), width = 0.7)+
    stat_summary(fun.data = mean_sdl,fun.args = list(mult=1),geom = "errorbar", width = 0.2)+
    stat_compare_means(comparisons = myComparision1, method = "t.test", label = "p.forma") +
    annotate("text", x=2, y= 0.89, label= round(sars_rawcounts_percentage_NT_mean,2)) + 
    annotate("text", x=3, y= 0.8, label= round(sars_rawcounts_percentage_T_mean,2)) +
    labs(x = "",y = "Covid19 Read Counts (%)")+
    theme_classic()
  p_covidReadCountsPercent
  
  #pdf("./Results/Figure/1_covidReadCounts.pdf", width = 4.0, height = 4.36)
  ggsave(filename = "./Results/Figure/2_covidReadCounts.pdf", p_covidReadCountsPercent, width = 4.0, height = 4.36)
  #dev.off()
}

# valcano of DE
{
  # DE Mock and NT
  data_DE_Mock_NT <- read.table("./Results/Table/DEG/DE_Data_Mock_NT.txt", header = T, sep = "\t", stringsAsFactors = F)
  # exclude sars genes
  data_DE_Mock_NT <- data_DE_Mock_NT[!data_DE_Mock_NT$Ensembl %in% sars_genes,]
  data_DE_Mock_NT$log2FoldChange <- -data_DE_Mock_NT$log2FoldChange
  p_valcano_Mock_NT <- qunplotValcano(dat = data_DE_Mock_NT, tagItem = "Mock VS NT")
  p_valcano_Mock_NT
  
  # DE NT and T
  data_DE_NT_T <- read.table("./Results/Table/DEG/DE_Data_NT_T.txt", header = T, sep = "\t", stringsAsFactors = F)
  # exclude sars genes
  data_DE_NT_T <- data_DE_NT_T[!data_DE_NT_T$Ensembl %in% sars_genes,]
  data_DE_NT_T$log2FoldChange <- -data_DE_NT_T$log2FoldChange
  p_valcano_NT_T <- qunplotValcano(dat = data_DE_NT_T, tagItem = "NT VS T")
  p_valcano_NT_T
}

# Expression Change
{
  # get up/down expressed genes of MockVSNT group
  Ensembl_Mock_NT_up <- setdiff(data_DE_Mock_NT[data_DE_Mock_NT$sig_strict == "up",]$Ensembl, sars_genes)
  Ensembl_Mock_NT_down <- setdiff(data_DE_Mock_NT[data_DE_Mock_NT$sig_strict == "down",]$Ensembl, sars_genes)

  exp_Mock_NT_up <- exp_data_fpkm[rownames(exp_data_fpkm) %in% Ensembl_Mock_NT_up,][,fpkm_mean_cols]
  exp_Mock_NT_up_scale <- data.frame(na.omit(t(scale(t(exp_Mock_NT_up),center = TRUE, scale = TRUE))))
  colnames(exp_Mock_NT_up_scale) <- c("Mock_mean_up", "NT_mean_up", "T_mean_up")
  
  exp_Mock_NT_down <- exp_data_fpkm[rownames(exp_data_fpkm) %in% Ensembl_Mock_NT_down,][,fpkm_mean_cols]
  exp_Mock_NT_down_scale <- data.frame(na.omit(t(scale(t(exp_Mock_NT_down), center = TRUE, scale = TRUE))))
  colnames(exp_Mock_NT_down_scale) <- c("Mock_mean_down", "NT_mean_down", "T_mean_down")
  
  pheatmap(na.omit(exp_Mock_NT_up_scale), cluster_cols = F, show_rownames = F, color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
  pheatmap(na.omit(exp_Mock_NT_down_scale), cluster_cols = F, show_rownames = F, color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
  
  exp_Mock_NT_up_scale$Check_Mock_NT_up <- 0
  exp_Mock_NT_up_scale$Check_Mock_T_up <- 0
  exp_Mock_NT_up_scale$Checkreverse <- 0
  
  for (i in 1:nrow(exp_Mock_NT_up_scale)){
    if (exp_Mock_NT_up_scale[i,]$NT_mean_up >= exp_Mock_NT_up_scale[i,]$Mock_mean_up){
      exp_Mock_NT_up_scale[i,]$Check_Mock_NT_up <- 1
    }
    if (exp_Mock_NT_up_scale[i,]$T_mean_up >= exp_Mock_NT_up_scale[i,]$Mock_mean_up){
      exp_Mock_NT_up_scale[i,]$Check_Mock_T_up <- 1
    }
    if(exp_Mock_NT_up_scale[i,]$T_mean_up <= exp_Mock_NT_up_scale[i,]$NT_mean_up){
      exp_Mock_NT_up_scale[i,]$Checkreverse <- 1
    }
  }
  
  exp_Mock_NT_down_scale$Check_Mock_NT_down <- 0
  exp_Mock_NT_down_scale$Check_Mock_T_down <- 0
  exp_Mock_NT_down_scale$Checkreverse <- 0
  
  for (i in 1:nrow(exp_Mock_NT_down_scale)){
    if (exp_Mock_NT_down_scale[i,]$NT_mean_down <= exp_Mock_NT_down_scale[i,]$Mock_mean_down){
      exp_Mock_NT_down_scale[i,]$Check_Mock_NT_down <- 1
    }
    if (exp_Mock_NT_down_scale[i,]$T_mean_down <= exp_Mock_NT_down_scale[i,]$Mock_mean_down){
      exp_Mock_NT_down_scale[i,]$Check_Mock_T_down <- 1
    }
    if(exp_Mock_NT_down_scale[i,]$T_mean_down >= exp_Mock_NT_down_scale[i,]$NT_mean_down){
      exp_Mock_NT_down_scale[i,]$Checkreverse <- 1
    }
  }
  
  # box plot
  exp_Mock_NT_up_scale_stack <- stack(exp_Mock_NT_up_scale[,1:3])
  exp_Mock_NT_down_scale_stack <- stack(exp_Mock_NT_down_scale[,1:3])
  exp_Mock_NT_scale_stack <- rbind(exp_Mock_NT_up_scale_stack, exp_Mock_NT_down_scale_stack)
  
  p_box_1 <- ggplot(exp_Mock_NT_up_scale_stack, aes(x=ind, y=values, fill = ind)) + 
    geom_boxplot(outlier.shape = NA) + 
    scale_fill_manual(values = c("grey", "Dark Red", "Brown")) +
    labs(x = "",y = "z-scaled fpkm")+
    theme_classic()
  p_box_1
  
  p_box_2 <- ggplot(exp_Mock_NT_down_scale_stack, aes(x=ind, y=values, fill = ind)) + 
    geom_boxplot(outlier.shape = NA) + 
    scale_fill_manual(values = c("grey","Dark Green", "Olive Drab")) +
    labs(x = "",y = "Covid19 Read Counts (%)")+
    theme_classic()
  p_box_2
  
  myComparision2 <- list(c("NT_mean_up","T_mean_up"), c("NT_mean_down","T_mean_down"))
  p_box_3 <- ggplot(exp_Mock_NT_scale_stack, aes(x=ind, y=values, fill = ind)) + 
    geom_boxplot(outlier.shape = NA) + 
    scale_fill_manual(values = c("grey", "Dark Red", "Brown", "black","Dark Green", "Olive Drab")) +
    #stat_summary(fun.data = mean_sdl,fun.args = list(mult=1), geom = "errorbar", width = 0.2)+
    stat_compare_means(comparisons = myComparision2, method = "wilcox.test", label = "p.forma") +
    labs(x = "",y = "Average FPKM (Z scaled)") +
    theme_classic() + 
    theme(axis.text.x = element_blank())
  p_box_3
  ggsave(filename = "./Results/Figure/3_expReverse.pdf", p_box_3, width = 5.6, height = 5.0)
}

# GO function analysis
{
  # Mock VS NT
  # dysregulated genes with pvalue <= 0.01, abs(log2FC) > 1
  
  
  # NT VS T
  
}


