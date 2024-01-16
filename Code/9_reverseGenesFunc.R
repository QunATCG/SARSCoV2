###
# @Descripttion: expression reverse genes func
# @Author: LiQun
# @Email: qun.li@ki.se
# @Date: 2 Jan 2024 (09:30:33)
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
  library(org.Hs.eg.db)
  library(reshape2)
  library(clusterProfiler)
  library(ggplot2)
  library(gplots)
  library(ggrepel)
  library(ggsci)
  library(ggpubr)
  library(UpSetR)
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
  
  # DEG table
  {
    deg_mock_nt_total <- read.table("./Results/Table/DEG/DE_Data_Mock_NT.txt", header = T, sep = "\t")
    deg_mock_t_total <- read.table("./Results/Table/DEG/DE_Data_Mock_T.txt", header = T, sep = "\t")
    deg_nt_t_total <- read.table("./Results/Table/DEG/DE_Data_NT_T.txt", header = T, sep = "\t")
  } 
    
  # load public exp
  {
    exp_GSE150316 <- read.table("./Results/Table/PublicExp/exp_GSE150316_scale.txt", header = T, sep = "\t")
    exp_GSE152418 <- read.table("./Results/Table/PublicExp/exp_GSE152418_scale.txt", header = T, sep = "\t")
    exp_GSE157103 <- read.table("./Results/Table/PublicExp/exp_GSE157103_scale.txt", header = T, sep = "\t")
  }
}

# reverse genes
{
  # creat reverse genes expression table
  NT_Mock_up_genes <- deg_mock_nt_total[deg_mock_nt_total$sig == "up",]$Ensembl
  NT_Mock_down_genes <- deg_mock_nt_total[deg_mock_nt_total$sig == "down",]$Ensembl
  NT_Mock_none_genes <- deg_mock_nt_total[deg_mock_nt_total$sig == "none",]$Ensembl
  
  T_Mock_up_genes <- deg_mock_t_total[deg_mock_t_total$sig == "up",]$Ensembl
  T_Mock_down_genes <- deg_mock_t_total[deg_mock_t_total$sig == "down",]$Ensembl
  T_Mock_none_genes <- deg_mock_t_total[deg_mock_t_total$sig == "none",]$Ensembl
  
  total_genes <- Reduce(union, list(NT_Mock_up_genes, NT_Mock_down_genes, T_Mock_up_genes, T_Mock_down_genes))
  both_up_genes <- intersect(NT_Mock_up_genes, T_Mock_up_genes)
  both_down_genes <- intersect(NT_Mock_down_genes, T_Mock_down_genes)
  
  total_genes_exp <- exp_fpkm[exp_fpkm$Ensembl %in% total_genes,]
  both_up_genes_exp <- exp_fpkm[exp_fpkm$Ensembl %in% both_up_genes,]
  both_down_genes_exp <- exp_fpkm[exp_fpkm$Ensembl %in% both_down_genes,]
  
  both_up_genes_ToverNT <- both_up_genes_exp[both_up_genes_exp$T_mean_scale >= both_up_genes_exp$NT_mean_scale, ]
  both_up_genes_TlessNT <- both_up_genes_exp[both_up_genes_exp$T_mean_scale < both_up_genes_exp$NT_mean_scale, ]
  
  both_down_genes_ToverNT <- both_down_genes_exp[both_down_genes_exp$T_mean_scale >= both_down_genes_exp$NT_mean_scale, ]
  both_down_genes_TlessNT <- both_down_genes_exp[both_down_genes_exp$T_mean_scale < both_down_genes_exp$NT_mean_scale, ]
  
  NT_unique_up_genes <- setdiff(NT_Mock_up_genes, T_Mock_up_genes)
  NT_unique_down_genes <- setdiff(NT_Mock_down_genes, T_Mock_down_genes)
  T_unique_up_genes <- setdiff(T_Mock_up_genes, NT_Mock_up_genes)
  T_unique_down_genes <- setdiff(T_Mock_down_genes, NT_Mock_down_genes)
  
  # output genes for DAVID
  
  markGene1 <- data.frame(Ensembl = union(both_up_genes_TlessNT$Ensembl, intersect(NT_Mock_up_genes, T_Mock_none_genes)))
  markGene1EG <- merge(markGene1, Ensembl_gene, by = "Ensembl")
  write.table(markGene1EG, "./Results/Table/DEG/DEG_both_up_genes_TlessNT.txt", row.names = F, quote = F, sep = "\t")
  
  markGene2 <- data.frame(Ensembl = union(both_down_genes_ToverNT$Ensembl, intersect(NT_Mock_down_genes, T_Mock_none_genes)))
  markGene2EG <- merge(markGene, Ensembl_gene, by = "Ensembl")
  
  write.table(markGeneEG, "./Results/Table/DEG/DEG_both_down_genes_ToverNT.txt", row.names = F, quote = F, sep = "\t")
  
  # GO analysis
  GO_NT_unique_up_genes <- as.data.frame(enrichGO(gene = NT_unique_up_genes, OrgDb = org.Hs.eg.db,
                                   keyType = 'ENSEMBL', ont = "BP", pAdjustMethod = "BH",
                                   pvalueCutoff = 0.05, qvalueCutoff = 0.1, readable = TRUE))
  GO_NT_unique_down_genes <- as.data.frame(enrichGO(gene = NT_unique_down_genes, OrgDb = org.Hs.eg.db,
                                                  keyType = 'ENSEMBL', ont = "BP", pAdjustMethod = "BH",
                                                  pvalueCutoff = 0.05, qvalueCutoff = 0.1, readable = TRUE))
  
  GO_T_unique_up_genes <- as.data.frame(enrichGO(gene = T_unique_up_genes, OrgDb = org.Hs.eg.db,
                                                  keyType = 'ENSEMBL', ont = "BP", pAdjustMethod = "BH",
                                                  pvalueCutoff = 0.05, qvalueCutoff = 0.1, readable = TRUE))
  GO_T_unique_down_genes <- as.data.frame(enrichGO(gene = T_unique_down_genes, OrgDb = org.Hs.eg.db,
                                                    keyType = 'ENSEMBL', ont = "BP", pAdjustMethod = "BH",
                                                    pvalueCutoff = 0.05, qvalueCutoff = 0.1, readable = TRUE))
  
  go_gene_1 <- markGene1EG[grep("^ENSG", markGene1EG$Gene, invert = T),]$Ensembl
  GO_both_up_genes_TlessNT <- as.data.frame(enrichGO(gene = go_gene_1, OrgDb = org.Hs.eg.db,
                                                     keyType = 'ENSEMBL', ont = "BP", pAdjustMethod = "BH",
                                                     pvalueCutoff = 0.05, qvalueCutoff = 0.1, readable = TRUE))
  
  go_gene_2 <- markGene2EG[grep("^ENSG", markGene2EG$Gene, invert = T),]$Ensembl
  GO_both_down_genes_ToverNT <- as.data.frame(enrichGO(gene = go_gene_2, OrgDb = org.Hs.eg.db,
                                                       keyType = 'ENSEMBL', ont = "BP", pAdjustMethod = "BH",
                                                       pvalueCutoff = 0.05, qvalueCutoff = 0.1, readable = TRUE))
  
  # creat upset input table
  upset_table_genes <- total_genes_exp[,c("Ensembl", "Gene", "Mock_mean_scale", "NT_mean_scale", "T_mean_scale")]
  upset_table_genes$NT_Mock_up <- 0
  upset_table_genes$NT_Mock_down <- 0
  upset_table_genes$T_Mock_up <- 0
  upset_table_genes$T_Mock_down <- 0
  for (i in 1:nrow(upset_table_genes)){
    if (upset_table_genes[i,]$Ensembl %in% NT_Mock_up_genes){
      upset_table_genes[i,]$NT_Mock_up <- 1
    }
    
    if (upset_table_genes[i,]$Ensembl %in% NT_Mock_down_genes){
      upset_table_genes[i,]$NT_Mock_down <- 1
    }
    
    if (upset_table_genes[i,]$Ensembl %in% T_Mock_up_genes){
      upset_table_genes[i,]$T_Mock_up <- 1
    }
    
    if (upset_table_genes[i,]$Ensembl %in% T_Mock_down_genes){
      upset_table_genes[i,]$T_Mock_down <- 1
    }
  }
  
  # def function for upset
  selectReverseGene <- function(row, geneList){
    data <- (row["Ensembl"] %in% geneList)
  }
  
  pdf("./Results/Figure/12_genesUpset.pdf", width = 6.4, height = 4.6)
  upset(upset_table_genes, nsets = 4, order.by = "freq", set_size.show = FALSE,
        main.bar.color = "grey", matrix.color = "black", sets.bar.color = "grey22",
        sets = c("NT_Mock_up", "NT_Mock_down", "T_Mock_up", "T_Mock_down"),
        queries = list(
          list(
            query = selectReverseGene,
            params = list(both_up_genes_TlessNT$Ensembl),
            color = "black",
            active = T
          ),
          list(
           query = selectReverseGene,
           params = list(both_down_genes_ToverNT$Ensembl),
           color = "black",
           active = T
          ),
          list(
          query = intersects, 
          params = list("NT_Mock_up", "T_Mock_up"), 
          color = "Dark Red", active = F
            ),
          list(
            query = intersects, 
            params = list("NT_Mock_down", "T_Mock_down"), 
            color = "Dark Green", active = F
          ))
        )
  dev.off()
}


# plot
{
  # detection of chemical stimulus involved in sensory perception of smell- GO:0050911
  genes_selected_0 <- c("OR1L4","OR5AU1","OR1L6","OR2AT4","OR10AD1","OR52W1","OR11H4","OR2B8P","OR1F2P","OR9A1P","OR9A4")
  genes_selected_fpkm <- exp_total[exp_total$Gene %in% genes_selected_0,][,c(11:20)]
  genes_selected_fpkm_melt <- melt(genes_selected_fpkm, var.id = Gene)

  mock_num <- length(grep("^Mock", genes_selected_fpkm_melt$variable))
  NT_num <- length(grep("^NT", genes_selected_fpkm_melt$variable))
  T_num <- length(grep("^T", genes_selected_fpkm_melt$variable))
  genes_selected_fpkm_melt$variable <- c(rep("Mock", mock_num), rep("NT", NT_num), rep("T", T_num))
  genes_selected_fpkm_melt$Gene <- factor(x= genes_selected_fpkm_melt$Gene, levels = c("OR9A1P","OR9A4","OR2B8P", "OR11H4","OR52W1", "OR1L4", "OR1L6", "OR5AU1","OR2AT4","OR10AD1","OR1F2P"))
  
  pdf("./Results/Figure/9_Histogram_mock_nt_up_smellGenes.pdf", width = 10.5, height = 3.7)
  ggplot(data = genes_selected_fpkm_melt, mapping = aes(x = Gene, y = value, fill = variable)) +
    geom_bar(stat = "summary", position = position_dodge(0.7)) +
    stat_summary(fun.data = "mean_sd", geom = "errorbar", colour = "black", width = 0.15, position = position_dodge(0.7)) +
    scale_fill_manual(values = c("#252525","#bdbdbd", "#525252")) +
    labs(x = "",y = "FPKM") +
    theme_classic()
  dev.off()

  # response to virus
  genes_selected_1 <- c("IFI6","RIGI","PLA2G10","IFIH1","IFIT3","IFIT2","OASL","IL23R","IFIT1", "MRC1")
  genes_selected_1_1 <- c("TPT1")
  genes_selected_fpkm <- exp_total[exp_total$Gene %in% genes_selected_1,][,c(11:20)]
  genes_selected_fpkm <- exp_total[exp_total$Gene %in% genes_selected_1_1,][,c(11:20)]
  genes_selected_fpkm_melt <- melt(genes_selected_fpkm, var.id = Gene)
  
  mock_num <- length(grep("^Mock", genes_selected_fpkm_melt$variable))
  NT_num <- length(grep("^NT", genes_selected_fpkm_melt$variable))
  T_num <- length(grep("^T", genes_selected_fpkm_melt$variable))
  genes_selected_fpkm_melt$variable <- c(rep("Mock", mock_num), rep("NT", NT_num), rep("T", T_num))
  genes_selected_fpkm_melt$Gene <- factor(x= genes_selected_fpkm_melt$Gene, levels = c("IFIT1","IFIT2", "IFIH1","IFIT3", "OASL", "RIGI", "IFI6", "IL23R","MRC1", "PLA2G10"))
  
  pdf("./Results/Figure/9_Histogram_nt_t_down_virusresponseGenes_1.pdf", width = 2.6, height = 3.7)
  ggplot(data = genes_selected_fpkm_melt, mapping = aes(x = Gene, y = value, fill = variable)) +
    geom_bar(stat = "summary", position = position_dodge(0.7)) +
    stat_summary(fun.data = "mean_sd", geom = "errorbar", colour = "black", width = 0.15, position = position_dodge(0.7)) +
    scale_fill_manual(values = c("#252525","#bdbdbd", "#525252")) +
    labs(x = "",y = "FPKM") +
    theme_classic()
  dev.off()
  
  
  
  
  # cytokines
  genes_seleted_cytokines_1 <- c("IL6", "TNF", "CCL7", "OSM", "TNFSF14", "GPLD1", "CLEC3B", "IFI44L", "OAS3")
  genes_seleted_cytokines_2 <- c("PRTN3", "LCN2", "CD24", "BPI", "DEFA4", "MMP8", "MPO", "TLR2")
  genes_seleted_cytokines_3 <- c("CISH", "TNFSF14", "OAS1", "GIMAP7", "GIMAP5", "TNFSF10", "DUSP2", "CXCR4", "RGS1", "TNFAIP3", "DUSP5", "AREG")
  
  
  
  # inflammatory response
  genes_selected_2 <- c("STAP1","PTGER3","IFI35","CXCL2","NLRC4","PTGS1","GGT1","NFKBIA","IL17F","HYAL1","SERPINC1",
                        "TNFAIP3","TEK","THEMIS2","SYT11","HRH4","NT5E","MYLK3","MEP1B","GPR32","ECM1","DUSP10",
                        "IL34","IL23R","CD200R1","PTX3","IL17RE","CXCL3","PTAFR","PTGER4","CSF1","LRRC19","CXCL17","LTA","GGT1")
  
  genes_selected_fpkm <- exp_fpkm[exp_fpkm$Gene %in% genes_selected_2,][,c(2,12:14)]
  rownames(genes_selected_fpkm) <- genes_selected_fpkm$Gene
  genes_selected_fpkm <- genes_selected_fpkm[,2:4]
  pheatmap(genes_selected_fpkm, scale = "row")
  
  # positive regulation of immune response
  genes_selected_3 <- c("CD74","BTN3A1","STAP1","IFI35","CEACAM1","NLRC4","TGFB2","LPXN","HPX","IL17F","PDE4D","HHLA2","COLEC11",
                        "THEMIS2","VTCN1","GPR32","CD226","IL23R","ERMAP","C2","PTAFR","GPR151","C1S","PLA2G6","FYB2","KLRK1","LTA","RAET1E")
  
  genes_selected_fpkm <- exp_fpkm[exp_fpkm$Gene %in% genes_selected_3,][,c(2,12:14)]
  rownames(genes_selected_fpkm) <- genes_selected_fpkm$Gene
  genes_selected_fpkm <- genes_selected_fpkm[,2:4]
  pheatmap(genes_selected_fpkm, scale = "row")
  
  # positive regulation of immune system process
  genes_selected_4 <- c("CD74","BTN3A1","STAP1","IFI35","CEACAM1","NLRC4","TGFB2","LPXN","HPX","IL17F","PDE4D","HES1","HHLA2","COLEC11",
                        "SOX4","THEMIS2","VTCN1","GPR32","DUSP10","CD226","JAM2","IL34","IL23R","ERMAP","C2","PTAFR","GPR151","JUN",
                        "HCLS1","C1S","CSF1","PLA2G6","FYB2","CXCL17","HLA-DMA","KLRK1","LTA","RAET1E")
  
  genes_selected_fpkm <- exp_fpkm[exp_fpkm$Gene %in% genes_selected_4,][,c(2,12:14)]
  rownames(genes_selected_fpkm) <- genes_selected_fpkm$Gene
  genes_selected_fpkm <- genes_selected_fpkm[,2:4]
  pheatmap(genes_selected_fpkm, scale = "row")
  
  # humoral immune response
  genes_selected_5 <- c("CXCL2","HPX","IL17F","COLEC11","IFNA5","CXCL3","C2","WFDC12","C1S","PLA2G6","LTA","H2BC6")
  genes_selected_fpkm <- exp_fpkm[exp_fpkm$Gene %in% genes_selected_5,][,c(2,12:14)]
  rownames(genes_selected_fpkm) <- genes_selected_fpkm$Gene
  genes_selected_fpkm <- genes_selected_fpkm[,2:4]
  pheatmap(genes_selected_fpkm, scale = "row")
  
  # genes imm
  genes_selected_imm <- Reduce(union, list(genes_selected_2, genes_selected_3, genes_selected_4, genes_selected_5))
  genes_selected_fpkm <- exp_fpkm[exp_fpkm$Gene %in% genes_selected_imm,][,c(2,12:14)]
  rownames(genes_selected_fpkm) <- genes_selected_fpkm$Gene
  genes_selected_fpkm <- genes_selected_fpkm[,2:4]
  colnames(genes_selected_fpkm) <- c("Mock", "NT", "T")
  
  p_pheatmap_immgenes_thisStudy <- pheatmap(genes_selected_fpkm, scale = "row", color = colorRampPalette(colors = c("#f0f0f0","#fff5f0","#67000d"))(100), 
           cellwidth = 10, cellheight = 5, cluster_cols = FALSE, clustering_method = "centroid", fontsize = 5, show_rownames = TRUE)
  
  ggsave(filename = "Results/Figure/10_mock_nt_up_immgenes_pheatmap.pdf", p_pheatmap_immgenes, width = 4.18, height = 13.55) 
  
  genes_selected_GSE157103 <- exp_GSE157103[exp_GSE157103$Gene %in% genes_selected_imm,][,c("Gene", "NC_Mean_scale", "Case_Mean_scale")]
  colnames(genes_selected_GSE157103) <- c("Gene", "NC3", "Case3")
  
  genes_selected_fpkm$Gene <- rownames(genes_selected_fpkm)
  gene_fpkm_GSE157103 <- merge(genes_selected_fpkm, genes_selected_GSE157103, by = "Gene")
  
  
  p_pheatmap_immgenes_thisStudy_GSE157103 <- pheatmap(gene_fpkm_GSE157103[,2:6], color = colorRampPalette(colors = c("#f0f0f0","#fff5f0","#67000d"))(100), 
                                            cellwidth = 10, cellheight = 2, cluster_cols = FALSE, clustering_method = "average", fontsize = 5, show_rownames = FALSE)
  
  ggsave(filename = "Results/Figure/10_mock_nt_up_immgenes_thisStudy_GSE157103_pheatmap.pdf", p_pheatmap_immgenes_thisStudy_GSE157103, width = 6.09, height = 3.44) 
  
  
  net_mock_nt_up_net  <- read.table("./Results/Table/Network/net_go_mock_NT_up_net.txt", header = T, sep = "\t")
  net_mock_nt_up_node <- read.table("./Results/Table/Network/net_go_mock_NT_up_node.txt", header = T, sep = "\t")
  p_mock_nt_up_immue <- qunplotmultiupbox(dat_net = net_mock_nt_up_node, exp1 = exp_fpkm, exp2 = exp_GSE152418, exp3 = exp_GSE157103, exp4 = exp_GSE150316, itemType = "immune")
  p_mock_nt_up_immue$p_thisStudy
  
  p_box_imm <- plot_grid(p_mock_nt_up_immue$p_thisStudy, p_mock_nt_up_immue$p_GSE157103, ncol = 2)
  ggsave(filename = "Results/Figure/10_mock_nt_up_immgenes_thisStudy_GSE157103_boxplot.pdf", p_box_imm, width = 6.09, height = 3.44) 
  
  targetGene <- "IFIT3"
  p_1 <- qunplotboxgenes(targetGene, exp_fpkm, myComparision = list(c("Mock", "NT"),c("Mock","T"), c("NT", "T")),myColors = c("grey", "Dark Red", "Brown"), isThisStudy = "T")
  p_1
  
  p_2 <- qunplotboxgenes(targetGene, exp_GSE150316, myComparision = list(c("NC", "Case")),myColors = c("grey", "Dark Red", "#fcbba1", "#fc9272"), isThisStudy = "316")
  p_2
  
  plot_grid(p_1, p_2, ncol = 2)
  
  p_3 <- qunplotboxgenes(targetGene, exp_GSE152418, myComparision = list(c("NC", "Case")),myColors = c("grey", "Dark Red", "Brown"), isThisStudy = "F")
  p_3
  
  p_4 <- qunplotboxgenes(targetGene, exp_GSE157103, myComparision = list(c("NC", "Case")),myColors = c("grey", "Dark Red", "Brown"), isThisStudy = "F")
  p_4

  # "BCL2A1" "CSF2" "EPSTI1" "MMP13" "CXCL6" "OAS2" "CXCL1" "CXCL2", "CXCL3", "IFI6", "IFI27" and "TNF"
  # (Gene expression profiling of SARS-CoV-2 infections reveal distinct primary lung cell and systemic immune infection responses that identify pathways relevant in COVID-19 disease)
  p_CCR_TNF <- qunplotboxgenes("TNF", exp_fpkm, myComparision = list(c("Mock", "NT"),c("Mock","T"), c("NT", "T")),myColors = c("grey", "Dark Red", "Brown"), isThisStudy = "T")
  
  # "CCL2" "CCL3" "CCL5" "CXCL8" "CXCL10"
  # (Chemokines and chemokine receptors during COVID-19 infection)
  p_CCR_CCR2 <- qunplotboxgenes("CCR2", exp_fpkm, myComparision = list(c("Mock", "NT"),c("Mock","T"), c("NT", "T")),myColors = c("grey", "Dark Red", "Brown"), isThisStudy = "T")
  p_CCR_CCR5 <- qunplotboxgenes("CCR5", exp_fpkm, myComparision = list(c("Mock", "NT"),c("Mock","T"), c("NT", "T")),myColors = c("grey", "Dark Red", "Brown"), isThisStudy = "T")
  
  p_CCR_CCL2 <- qunplotboxgenes("CCL2", exp_fpkm, myComparision = list(c("Mock", "NT"),c("Mock","T"), c("NT", "T")),myColors = c("grey", "Dark Red", "Brown"), isThisStudy = "T")
  p_CCR_CCL3 <- qunplotboxgenes("CCL3", exp_fpkm, myComparision = list(c("Mock", "NT"),c("Mock","T"), c("NT", "T")),myColors = c("grey", "Dark Red", "Brown"), isThisStudy = "T")
  
  p_CCR <- plot_grid(p_CCR_CCR2, p_CCR_CCR5, p_CCR_CCL2, p_CCR_CCL3, ncol = 4)
  ggsave(filename = "./Results/Figure/14_p_plotbox_CCR.pdf", p_CCR, width = 10, height = 3)
  
  # Cytokine pathway differences between low and high viral cases
  # Temporal and spatial heterogeneity of host response to SARS-CoV-2 pulmonary infection
  # JAK/STAT pathway
  p_1_JAK2 <- qunplotboxgenes("JAK2", exp_fpkm, myComparision = list(c("Mock", "NT"),c("Mock","T"), c("NT", "T")),myColors = c("grey", "Dark Red", "Brown"), isThisStudy = "T")
  p_2_JAK2 <- qunplotboxgenes("JAK2", exp_GSE150316, myComparision = list(c("NC", "Case")),myColors = c("grey", "Dark Red", "#fcbba1", "#fc9272"), isThisStudy = "316")
  
  p_1_STAT1 <- qunplotboxgenes("STAT1", exp_fpkm, myComparision = list(c("Mock", "NT"),c("Mock","T"), c("NT", "T")),myColors = c("grey", "Dark Red", "Brown"), isThisStudy = "T")
  p_2_STAT1 <- qunplotboxgenes("STAT1", exp_GSE150316, myComparision = list(c("NC", "Case")),myColors = c("grey", "Dark Red", "#fcbba1", "#fc9272"), isThisStudy = "316")
  
  #p_1_STAT2 <- qunplotboxgenes("STAT2", exp_fpkm, myComparision = list(c("Mock", "NT"),c("Mock","T"), c("NT", "T")),myColors = c("grey", "Dark Red", "Brown"), isThisStudy = "T")
  #p_2_STAT2 <- qunplotboxgenes("STAT2", exp_GSE150316, myComparision = list(c("NC", "Case")),myColors = c("grey", "Dark Red", "#fcbba1", "#fc9272"), isThisStudy = "316")
  #plot_grid(p_1_STAT2, p_2_STAT2)
  
  p_1_IL22RA2 <- qunplotboxgenes("IL22RA2", exp_fpkm, myComparision = list(c("Mock", "NT"),c("Mock","T"), c("NT", "T")),myColors = c("grey", "Dark Red", "Brown"), isThisStudy = "T")
  p_2_IL22RA2 <- qunplotboxgenes("IL22RA2", exp_GSE150316, myComparision = list(c("NC", "Case")),myColors = c("grey", "Dark Red", "#fcbba1", "#fc9272"), isThisStudy = "316")
  plot_grid(p_1_IL22RA2, p_2_IL22RA2)
  
  p_1_CXCR6 <- qunplotboxgenes("CXCR6", exp_fpkm, myComparision = list(c("Mock", "NT"),c("Mock","T"), c("NT", "T")),myColors = c("grey", "Dark Red", "Brown"), isThisStudy = "T")
  p_2_CXCR6 <- qunplotboxgenes("CXCR6", exp_GSE150316, myComparision = list(c("NC", "Case")),myColors = c("grey", "Dark Red", "#fcbba1", "#fc9272"), isThisStudy = "316")
  plot_grid(p_1_test, p_2_test)
  
  p_CP <- plot_grid(p_1_JAK2, p_2_JAK2, p_1_STAT1, p_2_STAT1,
            p_1_IL22RA2, p_2_IL22RA2, p_1_CXCR6, p_2_CXCR6,
            ncol = 4)
  
  ggsave(filename = "./Results/Figure/15_p_plotbox_CP.pdf", p_CP, width = 10, height = 6)
}



