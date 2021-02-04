library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)

batch1.markers <- FindAllMarkers(object = batch1, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
cluster_0 <- batch1.markers[batch1.markers$cluster == "0", ]
cluster_1 <- batch1.markers[batch1.markers$cluster == "1", ]
cluster_2 <- batch1.markers[batch1.markers$cluster == "2", ]
cluster_3 <- batch1.markers[batch1.markers$cluster == "3", ]
cluster_4 <- batch1.markers[batch1.markers$cluster == "4", ]
cluster_5 <- batch1.markers[batch1.markers$cluster == "5", ]
cluster_6 <- batch1.markers[batch1.markers$cluster == "6", ]
cluster_7 <- batch1.markers[batch1.markers$cluster == "7", ]
cluster_8 <- batch1.markers[batch1.markers$cluster == "8", ]

##### Cluster 0: Top 100 genes ######
cluster_0$Entrez.Gene <- mapIds(org.Hs.eg.db, keys=as.character(cluster_0$gene), 
                               keytype="SYMBOL", column="ENTREZID")

cluster_0 <- cluster_0[is.na(cluster_0$Entrez.Gene)==FALSE,]

cluster_0 <- cluster_0[!duplicated(cluster_0$Entrez.Gene),]

cluster_0 <- cluster_0[c(1:100),]

geneset <- as.character(cluster_0$Entrez.Gene)

# This will take a little while to run
ego <- enrichGO(gene = geneset,
                universe = NULL, #all available genes in database
                OrgDb = org.Hs.eg.db, #Hs: homo sapiens
                ont ="BP", #molecular function, biological process, cellular component
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,  # q value is FDR adjusted p value
                readable = TRUE) #will show gene symbol in images later rather than Entrez Gene ID
# dimensions - number of GO terms
dim(ego)

cluster0_GOgenes <- data.frame(ego$ID, ego$Description, ego$p.adjust, ego$geneID)

simp <- simplify(ego)

dim(simp)

cluster0_GOgenes <- data.frame(simp$ID, simp$Description, simp$p.adjust, simp$geneID)

barplot(simp)

emapplot(simp)

##### Cluster 1: Top 100 genes ######
cluster_1$Entrez.Gene <- mapIds(org.Hs.eg.db, keys=as.character(cluster_1$gene), 
                                keytype="SYMBOL", column="ENTREZID")

cluster_1 <- cluster_1[is.na(cluster_1$Entrez.Gene)==FALSE,]

cluster_1 <- cluster_1[!duplicated(cluster_1$Entrez.Gene),]

cluster_1 <- cluster_1[c(1:100),]

geneset <- as.character(cluster_1$Entrez.Gene)

# This will take a little while to run
ego <- enrichGO(gene = geneset,
                universe = NULL, #all available genes in database
                OrgDb = org.Hs.eg.db, #Hs: homo sapiens
                ont ="BP", #molecular function, biological process, cellular component
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,  # q value is FDR adjusted p value
                readable = TRUE) #will show gene symbol in images later rather than Entrez Gene ID
# dimensions - number of GO terms
dim(ego)

cluster1_GOgenes <- data.frame(ego$ID, ego$Description, ego$p.adjust, ego$geneID)

simp <- simplify(ego)

dim(simp)

cluster1_GOgenes <- data.frame(simp$ID, simp$Description, simp$p.adjust, simp$geneID)

barplot(simp)

emapplot(simp)

##### Cluster 2: Top 100 genes ######
cluster_2$Entrez.Gene <- mapIds(org.Hs.eg.db, keys=as.character(cluster_2$gene), 
                                keytype="SYMBOL", column="ENTREZID")

cluster_2 <- cluster_2[is.na(cluster_2$Entrez.Gene)==FALSE,]

cluster_2 <- cluster_2[!duplicated(cluster_2$Entrez.Gene),]

cluster_2 <- cluster_2[c(1:100),]

geneset <- as.character(cluster_2$Entrez.Gene)

# This will take a little while to run
ego <- enrichGO(gene = geneset,
                universe = NULL, #all available genes in database
                OrgDb = org.Hs.eg.db, #Hs: homo sapiens
                ont ="BP", #molecular function, biological process, cellular component
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,  # q value is FDR adjusted p value
                readable = TRUE) #will show gene symbol in images later rather than Entrez Gene ID
# dimensions - number of GO terms
dim(ego)

cluster2_GOgenes <- data.frame(ego$ID, ego$Description, ego$p.adjust, ego$geneID)

simp <- simplify(ego)

dim(simp)

cluster2_GOgenes <- data.frame(simp$ID, simp$Description, simp$p.adjust, simp$geneID)

barplot(simp)

emapplot(simp)

##### Cluster 3: Top 100 genes ######
cluster_3$Entrez.Gene <- mapIds(org.Hs.eg.db, keys=as.character(cluster_3$gene), 
                                keytype="SYMBOL", column="ENTREZID")

cluster_3 <- cluster_3[is.na(cluster_3$Entrez.Gene)==FALSE,]

cluster_3 <- cluster_3[!duplicated(cluster_3$Entrez.Gene),]

cluster_3 <- cluster_3[c(1:100),]

geneset <- as.character(cluster_3$Entrez.Gene)

# This will take a little while to run
ego <- enrichGO(gene = geneset,
                universe = NULL, #all available genes in database
                OrgDb = org.Hs.eg.db, #Hs: homo sapiens
                ont ="BP", #molecular function, biological process, cellular component
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,  # q value is FDR adjusted p value
                readable = TRUE) #will show gene symbol in images later rather than Entrez Gene ID
# dimensions - number of GO terms
dim(ego)

cluster3_GOgenes <- data.frame(ego$ID, ego$Description, ego$p.adjust, ego$geneID)

barplot(ego)


##### Cluster 4: Top 100 genes ######
cluster_4$Entrez.Gene <- mapIds(org.Hs.eg.db, keys=as.character(cluster_4$gene), 
                                keytype="SYMBOL", column="ENTREZID")

cluster_4 <- cluster_4[is.na(cluster_4$Entrez.Gene)==FALSE,]

cluster_4 <- cluster_4[!duplicated(cluster_4$Entrez.Gene),]

cluster_4 <- cluster_4[c(1:100),]

geneset <- as.character(cluster_4$Entrez.Gene)

# This will take a little while to run
ego <- enrichGO(gene = geneset,
                universe = NULL, #all available genes in database
                OrgDb = org.Hs.eg.db, #Hs: homo sapiens
                ont ="BP", #molecular function, biological process, cellular component
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,  # q value is FDR adjusted p value
                readable = TRUE) #will show gene symbol in images later rather than Entrez Gene ID
# dimensions - number of GO terms
dim(ego)

cluster4_GOgenes <- data.frame(ego$ID, ego$Description, ego$p.adjust, ego$geneID)

simp <- simplify(ego)

dim(simp)

cluster4_GOgenes <- data.frame(simp$ID, simp$Description, simp$p.adjust, simp$geneID)

barplot(simp)

emapplot(simp)

##### Cluster 5: Top 100 genes ######
cluster_5$Entrez.Gene <- mapIds(org.Hs.eg.db, keys=as.character(cluster_5$gene), 
                                keytype="SYMBOL", column="ENTREZID")

cluster_5 <- cluster_5[is.na(cluster_5$Entrez.Gene)==FALSE,]

cluster_5 <- cluster_5[!duplicated(cluster_5$Entrez.Gene),]

cluster_5 <- cluster_5[c(1:100),]

geneset <- as.character(cluster_5$Entrez.Gene)

# This will take a little while to run
ego <- enrichGO(gene = geneset,
                universe = NULL, #all available genes in database
                OrgDb = org.Hs.eg.db, #Hs: homo sapiens
                ont ="BP", #molecular function, biological process, cellular component
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,  # q value is FDR adjusted p value
                readable = TRUE) #will show gene symbol in images later rather than Entrez Gene ID
# dimensions - number of GO terms
dim(ego)

cluster5_GOgenes <- data.frame(ego$ID, ego$Description, ego$p.adjust, ego$geneID)

simp <- simplify(ego)

dim(simp)

cluster5_GOgenes <- data.frame(simp$ID, simp$Description, simp$p.adjust, simp$geneID)

barplot(simp)

emapplot(simp)

##### Cluster 6: Top 100 genes ######
cluster_6$Entrez.Gene <- mapIds(org.Hs.eg.db, keys=as.character(cluster_6$gene), 
                                keytype="SYMBOL", column="ENTREZID")

cluster_6 <- cluster_6[is.na(cluster_6$Entrez.Gene)==FALSE,]

cluster_6 <- cluster_6[!duplicated(cluster_6$Entrez.Gene),]

cluster_6 <- cluster_6[c(1:100),]

geneset <- as.character(cluster_6$Entrez.Gene)

# This will take a little while to run
ego <- enrichGO(gene = geneset,
                universe = NULL, #all available genes in database
                OrgDb = org.Hs.eg.db, #Hs: homo sapiens
                ont ="BP", #molecular function, biological process, cellular component
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,  # q value is FDR adjusted p value
                readable = TRUE) #will show gene symbol in images later rather than Entrez Gene ID
# dimensions - number of GO terms
dim(ego)

cluster6_GOgenes <- data.frame(ego$ID, ego$Description, ego$p.adjust, ego$geneID)

simp <- simplify(ego)

dim(simp)

cluster6_GOgenes <- data.frame(simp$ID, simp$Description, simp$p.adjust, simp$geneID)

barplot(simp)

emapplot(simp)

##### Cluster 7: Top 100 genes ######
cluster_7$Entrez.Gene <- mapIds(org.Hs.eg.db, keys=as.character(cluster_7$gene), 
                                keytype="SYMBOL", column="ENTREZID")

cluster_7 <- cluster_7[is.na(cluster_7$Entrez.Gene)==FALSE,]

cluster_7 <- cluster_7[!duplicated(cluster_7$Entrez.Gene),]

cluster_7 <- cluster_7[c(1:100),]

geneset <- as.character(cluster_7$Entrez.Gene)

# This will take a little while to run
ego <- enrichGO(gene = geneset,
                universe = NULL, #all available genes in database
                OrgDb = org.Hs.eg.db, #Hs: homo sapiens
                ont ="BP", #molecular function, biological process, cellular component
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,  # q value is FDR adjusted p value
                readable = TRUE) #will show gene symbol in images later rather than Entrez Gene ID
# dimensions - number of GO terms
dim(ego)

cluster7_GOgenes <- data.frame(ego$ID, ego$Description, ego$p.adjust, ego$geneID)

barplot(ego)

##### Cluster 8: Top 100 genes ######
cluster_8$Entrez.Gene <- mapIds(org.Hs.eg.db, keys=as.character(cluster_8$gene), 
                                keytype="SYMBOL", column="ENTREZID")

cluster_8 <- cluster_8[is.na(cluster_8$Entrez.Gene)==FALSE,]

cluster_8 <- cluster_8[!duplicated(cluster_8$Entrez.Gene),]

cluster_8 <- cluster_8[c(1:100),]

geneset <- as.character(cluster_8$Entrez.Gene)

# This will take a little while to run
ego <- enrichGO(gene = geneset,
                universe = NULL, #all available genes in database
                OrgDb = org.Hs.eg.db, #Hs: homo sapiens
                ont ="BP", #molecular function, biological process, cellular component
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,  # q value is FDR adjusted p value
                readable = TRUE) #will show gene symbol in images later rather than Entrez Gene ID
# dimensions - number of GO terms
dim(ego)

cluster8_GOgenes <- data.frame(ego$ID, ego$Description, ego$p.adjust, ego$geneID)

simp <- simplify(ego)

dim(simp)

cluster8_GOgenes <- data.frame(simp$ID, simp$Description, simp$p.adjust, simp$geneID)

barplot(simp)

emapplot(simp)

