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

# generation of network
{
  # select part of GO term to show
  go_mock_NT_up_show <- go_mock_NT_up[1:30,]
  net_go_mock_NT_up <- qunNetworkGeneration(go_mock_NT_up_show)
  
  go_mock_NT_down_show <- go_mock_NT_down[1:30,]
  net_go_mock_NT_down <- qunNetworkGeneration(go_mock_NT_down_show)

  net_go_NT_T_up <- qunNetworkGeneration(go_NT_T_up)
  net_go_NT_T_down <- qunNetworkGeneration(go_NT_T_down)
  
  # output for cytoscape
  write.table(net_go_mock_NT_down$net, "./Results/Table/Network/net_go_mock_NT_down_net.txt", quote = F, sep = "\t", row.names = F)
  #write.table(net_go_mock_NT_down$node, "./Results/Table/Network/net_go_mock_NT_down_node.txt", quote = F, sep = "\t", row.names = F)
  
  write.table(net_go_mock_NT_up$net, "./Results/Table/Network/net_go_mock_NT_up_net.txt", quote = F, sep = "\t", row.names = F)
  #write.table(net_go_mock_NT_up$node, "./Results/Table/Network/net_go_mock_NT_up_node.txt", quote = F, sep = "\t", row.names = F)
  
  write.table(net_go_NT_T_down$net, "./Results/Table/Network/net_go_NT_T_down_net.txt", quote = F, sep = "\t", row.names = F)
  #write.table(net_go_NT_T_down$node, "./Results/Table/Network/net_go_NT_T_down_node.txt", quote = F, sep = "\t", row.names = F)
  
  write.table(net_go_NT_T_up$net, "./Results/Table/Network/net_go_NT_T_up_net.txt", quote = F, sep = "\t", row.names = F)
  #write.table(net_go_NT_T_up$node, "./Results/Table/Network/net_go_NT_T_up_node.txt", quote = F, sep = "\t", row.names = F)
}