###
# @Descripttion: show gene expression
# @Author: LiQun
# @Email: qun.li@ki.se
# @Date: 30 Dec 2023 (14:27:27)
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

# libraries
{
  library(pheatmap)
  library(cowplot)
  library(ggplot2)
  library(gplots)
  library(ggrepel)
  library(ggsci)
  library(ggpubr)
  source("./Code/basicFunction.R")
}

# load data
{
  # Ensembl-gene pair
  Ensembl_gene <- read.table("./Data/SourceData/matchedID.txt", header = T, sep = "\t")
  
  # Expression
  # exp this study
  {
    exp_total <- read.table("./Data/ExpRNAseq/exp_Data_merge.txt", header = T, sep = "\t")
    exp_fpkm <- exp_total[,c(1,11:20)]
    exp_fpkm_scale <- t(scale(t(exp_fpkm[,c(3:11)])))
    exp_fpkm <- cbind(exp_fpkm[,c(1,2)], exp_fpkm_scale)
    exp_fpkm$Mock_mean_scale <- apply(exp_fpkm[,3:5], 1, mean)
    exp_fpkm$NT_mean_scale <- apply(exp_fpkm[,6:8], 1, mean)
    exp_fpkm$T_mean_scale <- apply(exp_fpkm[,9:11], 1, mean)
  }
  
  # exp public data
  # RNAseq analysis of PBMCs in a group of 17 COVID-19 subjects and 17 healthy controls
  {
    exp_GSE152418 <- read.table("./Data/SourceData/HumanPublicCovid19/GSE152418_p20047_Study1_RawCounts.txt", header = T, sep = "\t")
    colnames(exp_GSE152418)[1] <- "Ensembl"
    exp_GSE152418 <- merge(exp_GSE152418, Ensembl_gene, by = "Ensembl")
    exp_GSE152418 <- exp_GSE152418[!duplicated(exp_GSE152418$Gene),]
    exp_GSE152418 <- exp_GSE152418[,c(1,36,2:35)]
    colnames_exp_GSE152418 <- c("Ensembl", "Gene", paste("Case", 1:17, sep = "_"), paste("NC", 1:17, sep = "_"))
    colnames(exp_GSE152418) <- colnames_exp_GSE152418
    exp_GSE152418_scale <- t(scale(t(exp_GSE152418[,3:36])))
    exp_GSE152418 <- cbind(exp_GSE152418[,1:2], exp_GSE152418_scale)
    exp_GSE152418$Case_Mean_scale <- apply(exp_GSE152418[,c(3:19)],1,mean)
    exp_GSE152418$NC_Mean_scale <- apply(exp_GSE152418[,c(20:36)],1,mean)
    exp_GSE152418 <- na.omit(exp_GSE152418)
    write.table(exp_GSE152418, "./Results/Table/PublicExp/exp_GSE152418_scale.txt", sep = "\t", row.names = F, quote = F)
  }
  
  # 126 samples were analyzed in total with 100 COVID-19 patients and 26 non-COVID-19
  {
    exp_GSE157103 <- read.table("./Data/SourceData/HumanPublicCovid19/GSE157103_genes.tpm.tsv", header = T, sep = "\t")
    colnames(exp_GSE157103)[1] <- "Gene"
    exp_GSE157103 <- merge(exp_GSE157103, Ensembl_gene, by = "Gene")
    exp_GSE157103 <- exp_GSE157103[!duplicated(exp_GSE157103$Gene),]
    exp_GSE157103 <- exp_GSE157103[,c(128, 1, 2:127)]
    colnames_exp_GSE157103 <- c("Ensembl", "Gene", paste("Case", 1:100, sep = "_"), paste("NC", 1:26, sep = "_"))
    colnames(exp_GSE157103) <- colnames_exp_GSE157103
    exp_GSE157103_scale <- t(scale(t(exp_GSE157103[,3:128])))
    exp_GSE157103 <- cbind(exp_GSE157103[,c(1,2)], exp_GSE157103_scale)
    exp_GSE157103$Case_Mean_scale <- apply(exp_GSE157103[,c(3:102)], 1, mean)
    exp_GSE157103$NC_Mean_scale <- apply(exp_GSE157103[,c(103:128)], 1, mean)
    exp_GSE157103 <- na.omit(exp_GSE157103)
    write.table(exp_GSE157103, "./Results/Table/PublicExp/exp_GSE157103_scale.txt", sep = "\t", row.names = F, quote = F)
  }

  # # Autopsy samples from patients deceased due to SARS-Cov2 infection were collected for total RNA-seq analysis to assess viral load and immune response.
  # # lung/heart/liver/kidney/bowel/marrow
  {
    exp_GSE150316 <- read.table("./Data/SourceData/HumanPublicCovid19/GSE150316_RPMNormCounts_final.txt", header = T, sep = "\t")
    exp_GSE150316 <- merge(exp_GSE150316, Ensembl_gene, by = "Gene")
    exp_GSE150316 <- exp_GSE150316[!duplicated(exp_GSE150316$Gene),]
    exp_GSE150316 <- exp_GSE150316[,c(90,1,2:69)]
    colnames(exp_GSE150316) <- c("Ensembl","Gene","Case_1_High_lung1","Case_1_High_lung2","Case_1_High_lung3","Case_1_High_lung4","Case_1_High_heart1","Case_2_Low_lung1",
                                 "Case_2_Low_lung2","Case_2_Low_jejunum1","Case_2_Low_lung3","Case_2_Low_heart1","Case_3_Low_lung1","Case_3_Low_heart1","Case_3_Low_lung2",
                                 "Case_4_Low_lung1","Case_4_Low_heart1","Case_4_Low_liver1","Case_4_Low_lung2","Case_4_Low_kidney1","Case_4_Low_heart2","Case_4_Low_bowel1",
                                 "Case_5_High_lung1","Case_5_High_lung2","Case_5_High_lung3","Case_5_High_lung4","Case_5_High_fat1","Case_5_High_lung5","Case_5_High_skin1",
                                 "Case_5_High_bowel1","Case_5_High_liver1","Case_5_High_kidney1","Case_5_High_heart1","Case_5_High_marrow1","NC_1","NC_2",
                                 "NC_3","NC_4","NC_5","Case_6_Low_lung1","Case_6_Low_lung2","Case_6_Low_lung3","Case_6_Low_lung4","Case_6_Low_lung5",
                                 "Case_7_Low_lung1","Case_7_Low_lung2","Case_7_Low_lung3","Case_7_Low_lung4","Case_7_Low_lung5","Case_8_High_lung1","Case_8_High_lung2","Case_8_High_lung3",
                                 "Case_8_High_lung4","Case_8_High_lung5","Case_8_High_liver1","Case_8_High_bowel1","Case_8_High_heart1","Case_9_High_lung1","Case_9_High_lung2","Case_9_High_lung3",
                                 "Case_9_High_lung4","Case_9_High_lung5","Case_10_Low_lung1","Case_10_Low_lung2","Case_10_Low_lung3","Case_11_High_lung1","Case_11_High_lung2","Case_11_High_lung3",
                                 "Case_11_High_kidney1","Case_11_High_bowel1")
    exp_GSE150316 <- exp_GSE150316[,colnames(exp_GSE150316)[grep("Ensembl|Gene|NC|lung", colnames(exp_GSE150316))]]
    colcolnames(exp_GSE150316)[grep("High", colnames(exp_GSE150316))]
    exp_GSE150316_scale <- t(scale(t(exp_GSE150316[,c(3:49)])))
    exp_GSE150316 <- cbind(exp_GSE150316[,1:2], exp_GSE150316_scale)
    exp_GSE150316$Case_Mean_scale <- apply(exp_GSE150316[,c(3:18, 24:49)], 1, mean)
    exp_GSE150316$NC_Mean_scale <- apply(exp_GSE150316[,c(19:23)], 1, mean)
    exp_GSE150316$Case_High_Mean_scale <- apply(exp_GSE150316[,c(colnames(exp_GSE150316)[grep("High", colnames(exp_GSE150316))])], 1, mean)
    exp_GSE150316$Case_Low_Mean_scale <- apply(exp_GSE150316[,c(colnames(exp_GSE150316)[grep("Low", colnames(exp_GSE150316))])], 1, mean)
    exp_GSE150316 <- na.omit(exp_GSE150316)
    write.table(exp_GSE150316, "./Results/Table/PublicExp/exp_GSE150316_scale.txt", sep = "\t", row.names = F, quote = F)
  }
  
  
  # GO Network
  {
    # mock nt
    net_mock_nt_up_net  <- read.table("./Results/Table/Network/net_go_mock_NT_up_net.txt", header = T, sep = "\t")
    net_mock_nt_up_node <- read.table("./Results/Table/Network/net_go_mock_NT_up_node.txt", header = T, sep = "\t")
    
    net_mock_nt_down_net  <- read.table("./Results/Table/Network/net_go_mock_NT_down_net.txt", header = T, sep = "\t")
    net_mock_nt_down_node <- read.table("./Results/Table/Network/net_go_mock_NT_down_node.txt", header = T, sep = "\t")
    
    # nt t
    net_nt_t_up_net  <- read.table("./Results/Table/Network/net_go_NT_T_up_net.txt", header = T, sep = "\t")
    net_nt_t_up_node <- read.table("./Results/Table/Network/net_go_NT_T_up_node.txt", header = T, sep = "\t")
    
    net_nt_t_down_net  <- read.table("./Results/Table/Network/net_go_NT_T_down_net.txt", header = T, sep = "\t")
    net_nt_t_down_node <- read.table("./Results/Table/Network/net_go_NT_T_down_node.txt", header = T, sep = "\t")
  }
}

# plot
{
  ###############################
  # mock nt
  ###############################
  # mock nt up
  table(net_mock_nt_up_node$Type)
  # immue
  p_mock_nt_up_immue <- qunplotmultiupbox(dat_net = net_mock_nt_up_node, exp1 = exp_fpkm, exp2 = exp_GSE152418, exp3 = exp_GSE157103, exp4 = exp_GSE150316, itemType = "immune")
  
  # humoral
  p_mock_nt_up_humoral <- qunplotmultiupbox(dat_net = net_mock_nt_up_node, exp1 = exp_fpkm, exp2 = exp_GSE152418, exp3 = exp_GSE157103, exp4 = exp_GSE150316, itemType = "humoral")
  
  # inflammatory
  p_mock_nt_up_inflammatory <- qunplotmultiupbox(dat_net = net_mock_nt_up_node, exp1 = exp_fpkm, exp2 = exp_GSE152418, exp3 = exp_GSE157103, exp4 = exp_GSE150316, itemType = "inflammatory")
  
  # GPCR
  p_mock_nt_up_GPCR <- qunplotmultiupbox(dat_net = net_mock_nt_up_node, exp1 = exp_fpkm, exp2 = exp_GSE152418, exp3 = exp_GSE157103, exp4 = exp_GSE150316, itemType = "GPCR")
  
  # sensory perception
  p_mock_nt_up_SP <- qunplotmultiupbox(dat_net = net_mock_nt_up_node, exp1 = exp_fpkm, exp2 = exp_GSE152418, exp3 = exp_GSE157103, exp4 = exp_GSE150316, itemType = "sensory perception")
  
  # anion transport
  p_mock_nt_up_AT <- qunplotmultiupbox(dat_net = net_mock_nt_up_node, exp1 = exp_fpkm, exp2 = exp_GSE152418, exp3 = exp_GSE157103, exp4 = exp_GSE150316, itemType = "anion transport")
  
  # cell adhesion
  p_mock_nt_up_CA <- qunplotmultiupbox(dat_net = net_mock_nt_up_node, exp1 = exp_fpkm, exp2 = exp_GSE152418, exp3 = exp_GSE157103, exp4 = exp_GSE150316, itemType = "cell adhesion")
  
  # cilium movement
  p_mock_nt_up_CM <- qunplotmultiupbox(dat_net = net_mock_nt_up_node, exp1 = exp_fpkm, exp2 = exp_GSE152418, exp3 = exp_GSE157103, exp4 = exp_GSE150316, itemType = "cilium movement")
  
  # other
  p_mock_nt_up_other <- qunplotmultiupbox(dat_net = net_mock_nt_up_node, exp1 = exp_fpkm, exp2 = exp_GSE152418, exp3 = exp_GSE157103, exp4 = exp_GSE150316, itemType = "other")
  
  
  pdf("./Results/Figure/8_go_merge_mock_nt_up_public.pdf", width = 6.1, height = 14)
  plot_grid(p_mock_nt_up_immue$p_thisStudy, p_mock_nt_up_immue$p_GSE152418, p_mock_nt_up_immue$p_GSE157103, p_mock_nt_up_immue$p_GSE150316,
            p_mock_nt_up_humoral$p_thisStudy, p_mock_nt_up_humoral$p_GSE152418, p_mock_nt_up_humoral$p_GSE157103, p_mock_nt_up_humoral$p_GSE150316,
            p_mock_nt_up_inflammatory$p_thisStudy, p_mock_nt_up_inflammatory$p_GSE152418, p_mock_nt_up_inflammatory$p_GSE157103, p_mock_nt_up_inflammatory$p_GSE150316,
            p_mock_nt_up_GPCR$p_thisStudy, p_mock_nt_up_GPCR$p_GSE152418, p_mock_nt_up_GPCR$p_GSE157103, p_mock_nt_up_GPCR$p_GSE150316,
            p_mock_nt_up_SP$p_thisStudy, p_mock_nt_up_SP$p_GSE152418, p_mock_nt_up_SP$p_GSE157103, p_mock_nt_up_SP$p_GSE150316,
            p_mock_nt_up_AT$p_thisStudy, p_mock_nt_up_AT$p_GSE152418, p_mock_nt_up_AT$p_GSE157103, p_mock_nt_up_AT$p_GSE150316,
            p_mock_nt_up_CA$p_thisStudy, p_mock_nt_up_CA$p_GSE152418, p_mock_nt_up_CA$p_GSE157103, p_mock_nt_up_CA$p_GSE150316,
            p_mock_nt_up_CM$p_thisStudy, p_mock_nt_up_CM$p_GSE152418, p_mock_nt_up_CM$p_GSE157103, p_mock_nt_up_CM$p_GSE150316,
            p_mock_nt_up_other$p_thisStudy, p_mock_nt_up_other$p_GSE152418, p_mock_nt_up_other$p_GSE157103, p_mock_nt_up_other$p_GSE150316,
            ncol = 4)
  dev.off()
  # mock nt down
  table(net_mock_nt_down_node$Type)
  # apoptotic signaling pathway
  p_mock_nt_down_ASP <- qunplotmultidownbox(dat_net = net_mock_nt_down_node, exp1 = exp_fpkm, exp2 = exp_GSE152418, exp3 = exp_GSE157103, exp4 = exp_GSE150316, itemType = "apoptotic signaling pathway")
  
  # ATP
  p_mock_nt_down_ATP <- qunplotmultidownbox(dat_net = net_mock_nt_down_node, exp1 = exp_fpkm, exp2 = exp_GSE152418, exp3 = exp_GSE157103, exp4 = exp_GSE150316, itemType = "ATP")
  
  # Protein Fold
  p_mock_nt_down_PF <- qunplotmultidownbox(dat_net = net_mock_nt_down_node, exp1 = exp_fpkm, exp2 = exp_GSE152418, exp3 = exp_GSE157103, exp4 = exp_GSE150316, itemType = "Protein Fold")
  
  # regulation of proteolysis pathway
  p_mock_nt_down_RPP <- qunplotmultidownbox(dat_net = net_mock_nt_down_node, exp1 = exp_fpkm, exp2 = exp_GSE152418, exp3 = exp_GSE157103, exp4 = exp_GSE150316, itemType = "regulation of proteolysis pathway")
  
  # Translation
  p_mock_nt_down_Translation <- qunplotmultidownbox(dat_net = net_mock_nt_down_node, exp1 = exp_fpkm, exp2 = exp_GSE152418, exp3 = exp_GSE157103, exp4 = exp_GSE150316, itemType = "Translation")
  
  pdf("./Results/Figure/8_go_merge_mock_nt_down_public.pdf", width = 6.1, height = 5.5)
  plot_grid(p_mock_nt_down_ASP$p_thisStudy, p_mock_nt_down_ASP$p_GSE152418, p_mock_nt_down_ASP$p_GSE157103, p_mock_nt_down_ASP$p_GSE150316, 
            p_mock_nt_down_ATP$p_thisStudy, p_mock_nt_down_ATP$p_GSE152418, p_mock_nt_down_ATP$p_GSE157103,p_mock_nt_down_ATP$p_GSE150316,
            p_mock_nt_down_PF$p_thisStudy, p_mock_nt_down_PF$p_GSE152418, p_mock_nt_down_PF$p_GSE157103, p_mock_nt_down_PF$p_GSE150316,
            p_mock_nt_down_RPP$p_thisStudy, p_mock_nt_down_RPP$p_GSE152418, p_mock_nt_down_RPP$p_GSE157103, p_mock_nt_down_RPP$p_GSE150316,
            p_mock_nt_down_Translation$p_thisStudy, p_mock_nt_down_Translation$p_GSE152418, p_mock_nt_down_Translation$p_GSE157103, p_mock_nt_down_Translation$p_GSE150316,
            ncol = 4)
  dev.off()
  

  ###############################
  # nt t
  ###############################
  # nt t up
  table(net_nt_t_up_node$Type)
  # activation of adenylate cyclase activity
  p_nt_t_up_ACA <- qunplotmultiupbox(dat_net = net_nt_t_up_node, exp1 = exp_fpkm, exp2 = exp_GSE152418, exp3 = exp_GSE157103, exp4 = exp_GSE150316, itemType = "activation of adenylate cyclase activity")
  
  # chromatin assembly
  p_nt_t_up_CA <- qunplotmultiupbox(dat_net = net_nt_t_up_node, exp1 = exp_fpkm, exp2 = exp_GSE152418, exp3 = exp_GSE157103, exp4 = exp_GSE150316, itemType = "chromatin assembly")
  
  # splicing
  p_nt_t_up_splicing <- qunplotmultiupbox(dat_net = net_nt_t_up_node, exp1 = exp_fpkm, exp2 = exp_GSE152418, exp3 = exp_GSE157103, exp4 = exp_GSE150316, itemType = "Splicing")
  
  # synaptic signaling
  p_nt_t_up_SS <- qunplotmultiupbox(dat_net = net_nt_t_up_node, exp1 = exp_fpkm, exp2 = exp_GSE152418, exp3 = exp_GSE157103, exp4 = exp_GSE150316, itemType = "synaptic signaling")
  
  
  p_nt_t_up <- plot_grid(p_nt_t_up_ACA$p_thisStudy, p_nt_t_up_ACA$p_GSE152418, p_nt_t_up_ACA$p_GSE157103, p_nt_t_up_ACA$p_GSE150316,
                         p_nt_t_up_CA$p_thisStudy, p_nt_t_up_CA$p_GSE152418, p_nt_t_up_CA$p_GSE157103, p_nt_t_up_CA$p_GSE150316,
                         p_nt_t_up_splicing$p_thisStudy, p_nt_t_up_splicing$p_GSE152418, p_nt_t_up_splicing$p_GSE157103, p_nt_t_up_splicing$p_GSE150316,
                         p_nt_t_up_SS$p_thisStudy, p_nt_t_up_SS$p_GSE152418, p_nt_t_up_SS$p_GSE157103, p_nt_t_up_SS$p_GSE150316, ncol = 4)
  p_nt_t_up
  
  ggsave(filename = "./Results/Figure/8_go_merge_nt_t_up_public.pdf", p_nt_t_up, width = 11.34, height = 7.9)
  
  
  # nt t down
  table(net_nt_t_down_node$Type)
  # response to virus
  p_nt_t_down_virus <- qunplotmultidownbox(dat_net = net_nt_t_down_node, exp1 = exp_fpkm, exp2 = exp_GSE152418, exp3 = exp_GSE157103, exp4 = exp_GSE150316, itemType = "virus")
  
  p_nt_t_down_virus_merge <- plot_grid(p_nt_t_down_virus$p_thisStudy, p_nt_t_down_virus$p_GSE152418, p_nt_t_down_virus$p_GSE157103, p_nt_t_down_virus$p_GSE150316, ncol = 4)
  p_nt_t_down_virus_merge
  ggsave(filename = "./Results/Figure/8_go_merge_nt_t_down_public.pdf", p_nt_t_down_virus_merge, width = 11.34, height = 2)
  
  
  ###############################
  # nt t logFC
  ###############################
  # up
  DE_Data_NT_T <- read.table("./Results/Table/DEG/DE_Data_NT_T.txt", header = T, sep = "\t")
  # activation of adenylate cyclase activity
  dat_net_ACA_genes <- qunselectGenesBaseonGO(net_nt_t_up_node, "activation of adenylate cyclase activity")
  dat_net_ACA_genes_FC <- DE_Data_NT_T[DE_Data_NT_T$Gene %in% dat_net_ACA_genes,][,c("Gene", "log2FoldChange")]
  dat_net_ACA_genes_FC$Group <- "ACA"
  # synaptic signaling
  dat_net_SS_genes <- qunselectGenesBaseonGO(net_nt_t_up_node, "synaptic signaling")
  dat_net_SS_genes_FC <- DE_Data_NT_T[DE_Data_NT_T$Gene %in% dat_net_SS_genes,][,c("Gene", "log2FoldChange")]
  dat_net_SS_genes_FC$Group <- "SS"
  # chromatin assembly
  dat_net_CA_genes <- qunselectGenesBaseonGO(net_nt_t_up_node, "chromatin assembly")
  dat_net_CA_genes_FC <- DE_Data_NT_T[DE_Data_NT_T$Gene %in% dat_net_CA_genes,][,c("Gene", "log2FoldChange")]
  dat_net_CA_genes_FC$Group <- "CA"
  # Splicing
  dat_net_Splicing_genes <- qunselectGenesBaseonGO(net_nt_t_up_node, "Splicing")
  dat_net_Splicing_genes_FC <- DE_Data_NT_T[DE_Data_NT_T$Gene %in% dat_net_Splicing_genes,][,c("Gene", "log2FoldChange")]
  dat_net_Splicing_genes_FC$Group <- "Splicing"
  # down
  # response to virus
  dat_net_RV_genes <- qunselectGenesBaseonGO(net_nt_t_down_node, "virus")
  dat_net_RV_genes_FC <- DE_Data_NT_T[DE_Data_NT_T$Gene %in% dat_net_RV_genes,][,c("Gene", "log2FoldChange")]
  dat_net_RV_genes_FC$Group <- "virus"
  
  dat_net_selected_genes_FC <- Reduce(rbind, list(dat_net_ACA_genes_FC, dat_net_SS_genes_FC, dat_net_CA_genes_FC, dat_net_Splicing_genes_FC, dat_net_RV_genes_FC))
  dat_net_selected_genes_FC$Group <- factor(dat_net_selected_genes_FC$Group, levels = c("virus", "CA", "SS", "ACA","Splicing"))
  dat_net_selected_genes_FC$label = ifelse(dat_net_selected_genes_FC$Group == "virus", dat_net_selected_genes_FC$Gene, "")
  
  
  p_box_pathway <- ggplot(data = dat_net_selected_genes_FC, mapping = aes(x = Group, y = log2FoldChange, fill = Group)) +
    #stat_boxplot(geom = "errorbar", width = 0.15, aes(color = "black")) +
    geom_boxplot(size = 0.5, fill = "white", outlier.fill = "white", outlier.colour = "white") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    geom_jitter(aes(fill = Group), width = 0.2, shape = 21, size = 2.5) +
    scale_fill_manual(values = c("#991F1F", "black", "#28477D", "#50875C", "#CC4C02")) +
    ylab("log2 FoldChange (T / NT)") + xlab("") +
    geom_text_repel(data = dat_net_selected_genes_FC, aes(label = label), box.padding = unit(0.5, "lines"),
                    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE, size = 3) +
    theme_classic()
  
  ggsave(filename = "./Results/Figure/13_p_box_pathway.pdf", p_box_pathway, width = 5.4, height = 3.7)
}

