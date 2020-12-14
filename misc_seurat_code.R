#------------ Checking cluster genes from total data -----

setwd('/t1-data/user/lfelce/scRNA-Seq//t1-data/user/lfelce/scRNA-Seq/SmartSeq2_T-cells/')

tmarkers <- read.csv("Markers_for_Ling.csv", header=F)

cluster_0 <- tcell2.markers[tcell2.markers$cluster == "0", ]
cluster_1 <- tcell2.markers[tcell2.markers$cluster == "1", ]
cluster_2 <- tcell2.markers[tcell2.markers$cluster == "2", ]
cluster_3 <- tcell2.markers[tcell2.markers$cluster == "3", ]

cluster_0 <- cluster_0[cluster_0$p_val_adj < 0.05,]
cluster_1 <- cluster_1[cluster_1$p_val_adj < 0.05,]
cluster_2 <- cluster_2[cluster_2$p_val_adj < 0.05,]
cluster_3 <- cluster_3[cluster_3$p_val_adj < 0.05,]

genes_0 <- tmarkers[is.element(tmarkers$V1, cluster_0$gene),]
genes_1 <- tmarkers[is.element(tmarkers$V1, cluster_1$gene),]
genes_2 <- tmarkers[is.element(tmarkers$V1, cluster_2$gene),]
genes_3 <- tmarkers[is.element(tmarkers$V1, cluster_3$gene),]

markers_summary <- rbind(genes_0, genes_1, genes_2, genes_3)

#--------------- Checking cluster genes from CD4 data only -------

setwd('/t1-data/user/lfelce/scRNA-Seq//t1-data/user/lfelce/scRNA-Seq/SmartSeq2_T-cells/')

tmarkers <- read.csv("Markers_for_Ling.csv", header=F)

cluster_0 <- cd4.markers[cd4.markers$cluster == "0", ]
cluster_1 <- cd4.markers[cd4.markers$cluster == "1", ]
cluster_2 <- cd4.markers[cd4.markers$cluster == "2", ]

cluster_0 <- cluster_0[cluster_0$p_val_adj < 0.05,]
cluster_1 <- cluster_1[cluster_1$p_val_adj < 0.05,]
cluster_2 <- cluster_2[cluster_2$p_val_adj < 0.05,]

genes_0 <- tmarkers[is.element(tmarkers$V1, cluster_0$gene),]
genes_1 <- tmarkers[is.element(tmarkers$V1, cluster_1$gene),]
genes_2 <- tmarkers[is.element(tmarkers$V1, cluster_2$gene),]

markers_summary <- rbind(genes_0, genes_1, genes_2)


#########

test <- rep(c("A","B","C"), times =c("1","2","3"))

list <- c("CD4 CTL","CD4 Naive","CD4 Proliferating","CD4 TCM_1", "CD4_TCM_2", "CD4_TCM_3", "CD4_TEM_1","CD4_TEM2",
          "CD4_TEM_3", "CD4 TEM_4", "CD8 Naive", "CD8 Naive_2", "CD8 Proliferating", "CD8 TCM_1", "CD8 TCM_2",
          "CD8 TCM_3", "CD8 TEM_1","CD8 TEM_2","CD8 TEM_3","CD8 TEM_4","CD8 TEM_5","CD8 TEM_6", "MAIT", "Treg Memory",
          "Treg Naive")

cd4_ctl <- tmarkers[tmarkers$Cell.Type == "CD4 CTL", ]
cd4_naive <- tmarkers[tmarkers$Cell.Type == "CD4 Naive", ]
cd4_prof <- tmarkers[tmarkers$Cell.Type == "CD4 Proliferating", ]
cd4_tcm1 <- tmarkers[tmarkers$Cell.Type == "CD4 TCM_1", ]
cd4_tcm2 <- tmarkers[tmarkers$Cell.Type == "CD4 TCM_2", ]
cd4_tcm3 <- tmarkers[tmarkers$Cell.Type == "CD4 TCM_3", ]
cd4_tem1 <- tmarkers[tmarkers$Cell.Type == "CD4 TEM_1", ]
cd4_tem2 <- tmarkers[tmarkers$Cell.Type == "CD4 TEM_2", ]
cd4_tem3 <- tmarkers[tmarkers$Cell.Type == "CD4 TEM_3", ]
cd4_tem4 <- tmarkers[tmarkers$Cell.Type == "CD4 TEM_4", ]
cd8_naive <- tmarkers[tmarkers$Cell.Type == "CD8 Naive", ]
cd8_naive2 <- tmarkers[tmarkers$Cell.Type == "CD8 Naive_2", ]
cd8_prof <- tmarkers[tmarkers$Cell.Type == "CD8 Proliferating", ]
cd8_tcm1 <- tmarkers[tmarkers$Cell.Type == "CD8 TCM_1", ]
cd8_tcm2 <- tmarkers[tmarkers$Cell.Type == "CD8 TCM_2", ]
cd8_tcm3 <- tmarkers[tmarkers$Cell.Type == "CD8 TCM_3", ]
cd8_tem1 <- tmarkers[tmarkers$Cell.Type == "CD8 TEM_1", ]
cd8_tem2 <- tmarkers[tmarkers$Cell.Type == "CD8 TEM_2", ]
cd8_tem3 <- tmarkers[tmarkers$Cell.Type == "CD8 TEM_3", ]
cd8_tem4 <- tmarkers[tmarkers$Cell.Type == "CD8 TEM_4", ]
cd8_tem5 <- tmarkers[tmarkers$Cell.Type == "CD8 TEM_5", ]
cd8_tem6 <- tmarkers[tmarkers$Cell.Type == "CD8 TEM_6", ]
mait <- tmarkers[tmarkers$Cell.Type == "MAIT", ]
treg_mem <- tmarkers[tmarkers$Cell.Type == "Treg Memory", ]
treg_naive <- tmarkers[tmarkers$Cell.Type == "Treg Naive", ]
