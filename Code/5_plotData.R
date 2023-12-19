###
 # @Descripttion: plotData
 # @Author: Qun Li
 # @Email: qun.li@ki.se/liqun95@163.com
 # @Date: 2023-12-07 14:21:57
### 

# set working env
rm(list = ls())
# set your working path 
setwd("/Users/liqun/Desktop/Projects/Covid19/Data/Code/SARSCoV2/")
# check your working path
dirNow <- getwd()

### load libraries
library(ggplot2)
library(gplots)
library(ggsci)
library(ggpubr)
library(pheatmap)
library(VennDiagram)
library(clusterProfiler)
library(org.Hs.eg.db)


### function def
{
  keepOverExp <- function(dat, cutoff){
    res <- dat[which(rowMeans(dat[,c("Mock_1_FPKM", "Mock_2_FPKM", "Mock_3_FPKM")]) >=cutoff | 
                       rowMeans(dat[,c("T_1_FPKM", "T_2_FPKM", "T_3_FPKM")]) >=cutoff | 
                       rowMeans(dat[,c("NT_1_FPKM", "NT_2_FPKM", "NT_3_FPKM")]) >=cutoff),]
  }
  
  keepOverReadCount <- function(dat, cutoff, tagItem){
    if(tagItem == "star_total"){
      res <- dat[which(rowSums(dat[,c("Mock_1_ReadCount", "Mock_2_ReadCount", "Mock_3_ReadCount")] >= cutoff ) >= 3 |
                         rowSums(dat[,c("NT_1_ReadCount", "NT_2_ReadCount", "NT_3_ReadCount")] >= cutoff ) >= 3 |
                         rowSums(dat[,c("T_1_ReadCount", "T_2_ReadCount", "T_3_ReadCount")] >= cutoff ) >= 3),]
    }
    if(tagItem == "star_covid19"){
      res <- dat[which(rowSums(dat[,c("Mock_1_ExpectedCount", "Mock_2_ExpectedCount", "Mock_3_ExpectedCount")] >= cutoff ) >= 3 |
                         rowSums(dat[,c("NT_1_ExpectedCount", "NT_2_ExpectedCount", "NT_3_ExpectedCount")] >= cutoff ) >= 3 |
                         rowSums(dat[,c("T_1_ExpectedCount", "T_2_ExpectedCount", "T_3_ExpectedCount")] >= cutoff ) >= 3),]
    }
    if(tagItem == "hisat2_total"){
      res <- dat[which(rowSums(dat[,c("Mock_1_ReadCount", "Mock_2_ReadCount", "Mock_3_ReadCount")] >= cutoff ) >= 3 |
                         rowSums(dat[,c("NT_1_ReadCount", "NT_2_ReadCount", "NT_3_ReadCount")] >= cutoff ) >= 3 |
                         rowSums(dat[,c("T_1_ReadCount", "T_2_ReadCount", "T_3_ReadCount")] >= cutoff ) >= 3),]
    }
    if(tagItem == "hisat2_covid19"){
      res <- dat[which(rowSums(dat[,c("Mock_1_ReadCount", "Mock_2_ReadCount", "Mock_3_ReadCount")] >= cutoff ) >= 3 |
                         rowSums(dat[,c("NT_1_ReadCount", "NT_2_ReadCount", "NT_3_ReadCount")] >= cutoff ) >= 3 |
                         rowSums(dat[,c("T_1_ReadCount", "T_2_ReadCount", "T_3_ReadCount")] >= cutoff ) >= 3),]
    }
    return(res)
  }
  
  plotValcano <- function(dat, tagItem, baseline){
    dat <- na.omit(dat)
    ggplot(dat,aes(log2FoldChange,-log10(padj),color = sig))+ 
      geom_point()+
      scale_color_manual(values = c(down = "#00A087B2", up = "#DC0000B2", none = "grey")) +
      labs(x=expression(Log[2]*" Fold Change"), y=expression(-Log[10]*" (padj)")) +
      annotate("text", x=0, y= (baseline + 20), label= tagItem) + 
      annotate("text", x=0, y= (baseline + 10), label= paste("up:",table(dat$sig)[3])) + 
      annotate("text", x=0, y= baseline, label= paste("down:",table(dat$sig)[1])) +
      theme_classic()
  }
  
  qunGO <- function(datGene){
    GeneForGO <- datGene
    ego_ALL <- enrichGO(gene = GeneForGO, 
                        #universe = dataForGO$Ensembl,
                        OrgDb = org.Hs.eg.db,
                        keyType = 'ENSEMBL',
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05,
                        readable = TRUE)
    res <- as.data.frame(ego_ALL)
    return(res)
  }
}

### load data
### calculate covid19 readcounts
### hisat2
{
  countsData_total <- read.table("./Data/ExpRNAseq/Hisat2/Hisat2_exp_total.txt", header = T, sep = "\t")[,c("Ensembl", "Gene",
                      "Mock_1_ReadCount", "Mock_2_ReadCount", "Mock_3_ReadCount", "NT_1_ReadCount", "NT_2_ReadCount",
                      "NT_3_ReadCount", "T_1_ReadCount", "T_2_ReadCount", "T_3_ReadCount")]
  countsData_covid <- read.table("./Data/ExpRNAseq/Hisat2/Hisat2_exp_Covid19.txt", header = T, sep = "\t")[,c("Ensembl", "Gene",
                       "Mock_1_ReadCount", "Mock_2_ReadCount", "Mock_3_ReadCount", "NT_1_ReadCount", "NT_2_ReadCount",
                       "NT_3_ReadCount", "T_1_ReadCount", "T_2_ReadCount", "T_3_ReadCount")]
  Covid19_geneName <- c("ORF1ab","ORF1a","S", "ORF3a", "E", "M", "ORF6", "ORF7a","ORF7b",
                        "ORF8", "N", "ORF10")
  Covid19_Ensembl <- countsData_covid$Ensembl
  countsData_covid$Gene <- Covid19_geneName
  countsData <- rbind(countsData_total, countsData_covid)
  totalReadCounts <- colSums(countsData[,c("Mock_1_ReadCount", "Mock_2_ReadCount", "Mock_3_ReadCount", "NT_1_ReadCount", "NT_2_ReadCount",
                                           "NT_3_ReadCount", "T_1_ReadCount", "T_2_ReadCount", "T_3_ReadCount")])
  covidReadCounts <- colSums(countsData[c(57050:57061),c("Mock_1_ReadCount", "Mock_2_ReadCount", "Mock_3_ReadCount", "NT_1_ReadCount", "NT_2_ReadCount",
                                                         "NT_3_ReadCount", "T_1_ReadCount", "T_2_ReadCount", "T_3_ReadCount")])
  covidReadCountsPercent <- covidReadCounts/totalReadCounts
  covidReadCountsPercentDataFrame1 <- data.frame(sample = rep(c("Mock", "NT", "T"), each = 3), value = covidReadCountsPercent)
}


### star
#{
#  countsData_total <- read.table("./Data/ExpRNAseq/STAR/Star_exp_total.txt", header = T, sep = "\t")[,c("Ensembl", "Gene",
#                                                                                                        "Mock_1_ExpectedCount", "Mock_2_ExpectedCount", "Mock_3_ExpectedCount", "NT_1_ExpectedCount", "NT_2_ExpectedCount",
#                                                                                                        "NT_3_ExpectedCount", "T_1_ExpectedCount", "T_2_ExpectedCount", "T_3_ExpectedCount")]
#  countsData_covid <- read.table("./Data/ExpRNAseq/STAR/Star_exp_Covid19.txt", header = T, sep = "\t")[,c("Ensembl", "Gene",
#                                                                                                          "Mock_1_ExpectedCount", "Mock_2_ExpectedCount", "Mock_3_ExpectedCount", "NT_1_ExpectedCount", "NT_2_ExpectedCount",
#                                                                                                          "NT_3_ExpectedCount", "T_1_ExpectedCount", "T_2_ExpectedCount", "T_3_ExpectedCount")]
#  countsData <- rbind(countsData_total, countsData_covid)
#  totalReadCounts <- colSums(countsData[,c("Mock_1_ExpectedCount", "Mock_2_ExpectedCount", "Mock_3_ExpectedCount", "NT_1_ExpectedCount", "NT_2_ExpectedCount",
#                                           "NT_3_ExpectedCount", "T_1_ExpectedCount", "T_2_ExpectedCount", "T_3_ExpectedCount")])
#  covidReadCounts <- colSums(countsData[c(62755:62766),c("Mock_1_ExpectedCount", "Mock_2_ExpectedCount", "Mock_3_ExpectedCount", "NT_1_ExpectedCount", "NT_2_ExpectedCount",
#                                                         "NT_3_ExpectedCount", "T_1_ExpectedCount", "T_2_ExpectedCount", "T_3_ExpectedCount")])
#  covidReadCountsPercent <- covidReadCounts/totalReadCounts
#  covidReadCountsPercentDataFrame2 <- data.frame(sample = rep(c("Mock", "NT", "T"), each = 3), value = covidReadCountsPercent)
#}

NTMean <- mean(covidReadCountsPercentDataFrame1$value[4:6])
TMean <- mean(covidReadCountsPercentDataFrame1$value[7:9])

pvalueNT_T <- t.test(covidReadCountsPercentDataFrame1$value[4:6], covidReadCountsPercentDataFrame1$value[7:9])$p.value
myComparision <- list(c("NT","T"))

pbar_covidReadCountsPercent <- ggplot(covidReadCountsPercentDataFrame1,aes(sample,value))  +
  stat_summary(mapping=aes(fill = sample),fun=mean,geom = "bar",fun.args = list(mult=1),width=0.7)+
  stat_summary(fun.data=mean_sdl,fun.args = list(mult=1),geom="errorbar",width=0.2)+
  stat_compare_means(comparisons = myComparision, method = "t.test", label = "p.forma") +
  annotate("text", x=2, y= 0.89, label= round(NTMean,2)) + 
  annotate("text", x=3, y= 0.8, label= round(TMean,2)) +
  labs(x = "",y = "Covid19 Read Counts (%)")+
  #scale_y_continuous(expand = c(0,0),limits = c(0,2))+
  theme_classic() #+
  #theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
  #      axis.line=element_line(colour="black",size=0.25),
  #      axis.title=element_text(size=13,color="black"),
  #      axis.text = element_text(size=12,color="black"),
  #      legend.position="none")

pdf("./Results/Figure/1_covidReadCounts.pdf", width = 4.0, height = 4.36)
pbar_covidReadCountsPercent
dev.off()


### Only read hisat2
exp_NT_Mock <- keepOverExp(keepOverReadCount(read.csv("./Results/Table/DEG/Hisat2/DEG_NT_Mock_hisat2_result.csv", 
                               header = T, sep = ";", stringsAsFactors = FALSE, dec = ","), 5, "hisat2_total"),1)
#exp_NT_Mock <- subset(exp_NT_Mock, sig != "none")
go_NT_Mock_down <- read.csv("./Results/Table/GO/Hisat2/DEG_NT_Mock_hisat2_GO_Down.csv", header = T, sep = ";")
go_NT_Mock_up <- read.csv("./Results/Table/GO/Hisat2/DEG_NT_Mock_hisat2_GO_Up.csv", header = T, sep = ";")


exp_NT_T <- keepOverExp(keepOverReadCount(read.csv("./Results/Table/DEG/Hisat2/DEG_NT_T_hisat2_result.csv", 
                               header = T, sep = ";", stringsAsFactors = FALSE, dec = ","), 5, "hisat2_total"),1)
#exp_NT_T <- subset(exp_NT_T, sig != "none")
go_NT_T_down <- read.csv("./Results/Table/GO/Hisat2/DEG_NT_T_hisat2_GO_Down.csv", header = T, sep = ";")
go_NT_T_up <- read.csv("./Results/Table/GO/Hisat2/DEG_NT_T_hisat2_GO_Up.csv", header = T, sep = ";")


exp_T_Mock <- keepOverExp(keepOverReadCount(read.csv("./Results/Table/DEG/Hisat2/DEG_T_Mock_hisat2_result.csv",
                                header = T, sep = ";", stringsAsFactors = FALSE, dec = ","), 5, "hisat2_total"), 1)
#exp_T_Mock <- subset(exp_T_Mock, sig != "none")
go_T_Mock_down <- read.csv("./Results/Table/GO/Hisat2/DEG_T_Mock_hisat2_GO_Down.csv", header = T, sep = ";")
go_T_Mock_up <- read.csv("./Results/Table/GO/Hisat2/DEG_T_Mock_hisat2_GO_Up.csv", header = T, sep = ";")

### show DEG genes
### NT VS Mock
### exclude Covid genes
p_DEG_Valcano_NT_Mock <- plotValcano(exp_NT_Mock[13:nrow(exp_NT_Mock),],"NT VS Mock", 120)
pdf("./Results/Figure/2_p_DEG_Valcano_NT_Mock.pdf", width = 5.8, height = 4.6)
p_DEG_Valcano_NT_Mock
dev.off()

p_DEG_Valcano_T_Mock <- plotValcano(exp_T_Mock[13:nrow(exp_T_Mock),],"T VS Mock", 120)
pdf("./Results/Figure/2_p_DEG_Valcano_T_Mock.pdf", width = 5.8, height = 4.6)
p_DEG_Valcano_T_Mock
dev.off()

p_DEG_Valcano_T_NT <- plotValcano(exp_NT_T[13:nrow(exp_NT_Mock),],"T VS NT", 30)
pdf("./Results/Figure/2_p_DEG_Valcano_T_NT.pdf", width = 5.8, height = 4.6)
p_DEG_Valcano_T_NT
dev.off()

NT_Mock_up_genes <- na.omit(setdiff(exp_NT_Mock[exp_NT_Mock$sig == "up",]$Ensembl, Covid19_Ensembl))
NT_Mock_down_genes <- na.omit(setdiff(exp_NT_Mock[exp_NT_Mock$sig == "down",]$Ensembl, Covid19_Ensembl))

T_Mock_up_genes <- na.omit(setdiff(exp_T_Mock[exp_T_Mock$sig == "up",]$Ensembl, Covid19_Ensembl))
T_Mock_down_genes <- na.omit(setdiff(exp_T_Mock[exp_T_Mock$sig == "down",]$Ensembl, Covid19_Ensembl))

T_NT_up_genes <- na.omit(setdiff(exp_NT_T[exp_NT_T$sig == "up",]$Ensembl, Covid19_Ensembl))
T_NT_down_genes <- na.omit(setdiff(exp_NT_T[exp_NT_T$sig == "down",]$Ensembl, Covid19_Ensembl))

write.table(NT_Mock_up_genes, "./Results/Table/NT_Mock_up_genes.txt", quote = F, row.names = F, col.names = F)
write.table(NT_Mock_down_genes, "./Results/Table/NT_Mock_down_genes.txt", quote = F, row.names = F, col.names = F)
write.table(T_Mock_up_genes, "./Results/Table/T_Mock_up_genes.txt", quote = F, row.names = F, col.names = F)
write.table(T_Mock_down_genes, "./Results/Table/T_Mock_down_genes.txt", quote = F, row.names = F, col.names = F)
write.table(T_NT_up_genes, "./Results/Table/T_NT_up_genes.txt", quote = F, row.names = F, col.names = F)
write.table(T_NT_down_genes, "./Results/Table/T_NT_down_genes.txt", quote = F, row.names = F, col.names = F)

venn.diagram(
  x = list(NT_Mock_up_genes, NT_Mock_down_genes, T_Mock_up_genes, T_Mock_down_genes),
  category.names = c("NT_Mock_up" , "NT_Mock_down" , "T_Mock_up", "T_Mock_down"),
  filename = './Results/Figure/3_DEG_Overlap_1.tiff', imagetype = "tiff",
  output=T, disable.logging = T,  width = 4000
)

venn.diagram(
  x = list(NT_Mock_up_genes, NT_Mock_down_genes, T_NT_up_genes, T_NT_down_genes),
  category.names = c("NT_Mock_up" , "NT_Mock_down" , "T_NT_up", "T_NT_down"),
  filename = './Results/Figure/3_DEG_Overlap_2.tiff', imagetype = "tiff",
  output=T, disable.logging = T,  width = 4000
)


### merge expdata
exp_Data <- rbind(subset(exp_NT_Mock[,1:36], exp_NT_Mock[,1:36]$sig != "none"), 
                  subset(exp_NT_T[,1:36], exp_NT_T[,1:36]$sig != "none"), 
                  subset(exp_T_Mock[,1:36], exp_T_Mock[,1:36]$sig != "none"))
exp_Data <- exp_Data[!duplicated(exp_Data$Ensembl,),]
exp_Data <- exp_Data[!duplicated(exp_Data$Gene,),]

T_Uniq_up_genes <- setdiff(T_NT_up_genes, NT_Mock_up_genes)
T_Uniq_down_genes <- setdiff(T_NT_down_genes, NT_Mock_down_genes)
T_Uniq_genes <- c(T_Uniq_up_genes, T_Uniq_down_genes)
T_Uniq_genes_exp <- exp_Data[exp_Data$Ensembl %in% T_Uniq_genes,]

T_Uniq_genes_exp_Plot <- T_Uniq_genes_exp[,8:16]
rownames(T_Uniq_genes_exp_Plot) <- T_Uniq_genes_exp$Gene

pdf("./Results/Figure/4_DEG_Uniq_T.pdf", width = 5.69, height = 11.9)
pheatmap(T_Uniq_genes_exp_Plot, scale = "row", show_rownames = T, cluster_cols = F, cutree_rows = 2)
dev.off()

kk_total <- qunGO(T_Uniq_genes_exp$Ensembl)
kk_up <- qunGO(T_Uniq_up_genes)
kk_down <- qunGO(T_Uniq_down_genes)



