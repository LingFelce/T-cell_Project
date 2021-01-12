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

genes_0 <- tmarkers[is.element(tmarkers$V1, cluster_0$gene),]
genes_1 <- tmarkers[is.element(tmarkers$V1, cluster_1$gene),]
genes_2 <- tmarkers[is.element(tmarkers$V1, cluster_2$gene),]
genes_3 <- tmarkers[is.element(tmarkers$V1, cluster_3$gene),]
genes_4 <- tmarkers[is.element(tmarkers$V1, cluster_4$gene),]

markers_summary <- rbind(genes_0, genes_1, genes_2)


#----------------------- Separate patients and then integrate properly----------------
# https://satijalab.org/seurat/v3.2/integration.html
# split by patient - orig.ident from metadata
patient.list <- SplitObject(tcell2, split.by="orig.ident")
patient.list <- patient.list[c("005", "022", "025", "1062", "1105", "1131-TP-1", "1131-TP-2", 
                               "1134-TP-2", "1153", "1201-TP-2", "1493", "1504", "1525",
                               "1525-TP-1", "1525-TP-2")]

for (i in 1:length(patient.list)) {
  patient.list[[i]] <- NormalizeData(patient.list[[i]], verbose = FALSE)
  patient.list[[i]] <- FindVariableFeatures(patient.list[[i]], mean.function = ExpMean, 
                                            dispersion.function = LogVMR, 
                                            x.low.cutoff = 0.0125, x.high.cutoff = 3, 
                                            y.cutoff = 0.5, nfeatures = 2000)
}


# dimensions 10 to 50 is ok
reference.list <- patient.list[c("022", "025", "1062", "1493", "1504", "1525")]
patient.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)

patient.integrated <- IntegrateData(anchorset = patient.anchors, dims = 1:30)


#--------------------- CD8 separate patients and integrate properly ----------------
cd8_fresh <- CreateSeuratObject(counts = cd8_data, 
                                min.cells = 3, min.features = 200, 
                                project = "CD8 T-cell_data", assay = "RNA")

# split by patient - orig.ident from metadata (new patient 1131-TP-1)
cd8.list <- SplitObject(cd8_fresh, split.by="orig.ident")
cd8.list <- cd8.list[c("005", "1105", "1131-TP-1", "1131-TP-2", 
                               "1134-TP-2", "1153", "1201-TP-2",
                               "1525-TP-1", "1525-TP-2")]

for (i in 1:length(cd8.list)) {
  cd8.list[[i]] <- SCTransform(cd8.list[[i]], verbose = FALSE)
}

cd8.features <- SelectIntegrationFeatures(object.list = cd8.list, nfeatures = 2000)

cd8.list <- PrepSCTIntegration(object.list = cd8.list, anchor.features = cd8.features, 
                                    verbose = FALSE)


cd8.anchors <- FindIntegrationAnchors(object.list = cd8.list, normalization.method = "SCT", 
                                           anchor.features = cd8.features, 
                                      k.anchor = 5,
                                      k.filter = 46,
                                      k.score = 30)

cd8.integrated <- IntegrateData(anchorset = cd8.anchors, normalization.method = "SCT", k.weight=46)


cd8.integrated <- RunPCA(cd8.integrated, npcs = 30, verbose = FALSE)

cd8.integrated <- RunUMAP(cd8.integrated, reduction = "pca", dims = 1:10, verbose = FALSE)


Idents(cd8.integrated) <- "seurat_clusters"
UMAP_int_cluster <- DimPlot(cd8.integrated, reduction = "umap", group.by = "seurat_clusters")

UMAP_int_patient <- DimPlot(cd8.integrated, reduction = "umap", group.by = "orig.ident")

UMAP_int_cluster + UMAP_int_patient


#########
# Data wrangling to get bulk samples into 1 counts table/expression matrix

setwd('/t1-data/user/lfelce/scRNA-Seq/SmartSeq2_T-cells/')

## load data (from featureCounts)
all_201203 <-  fread('201203_counts.txt', stringsAsFactors = F, header=T)

# remove columns with chromosome, start, end, strand and length info
all_201203 <- all_201203[,-c(2:6)]

# make Geneid into row names
all_201203 <-tibble::column_to_rownames(all_201203, "Geneid")

# tidy up sample names
names(all_201203) <- gsub(x = names(all_201203), pattern = "./", replacement = "")  
names(all_201203) <- gsub(x = names(all_201203), pattern = ".bam", replacement = "")

# list of sample names
list_201203 <- as.data.frame(colnames(all_201203))

## load data (from featureCounts)
all_201221 <-  fread('201221_counts.txt', stringsAsFactors = F, header=T)

# remove columns with chromosome, start, end, strand and length info
all_201221 <- all_201221[,-c(2:6)]

# make Geneid into row names
all_201221 <-tibble::column_to_rownames(all_201221, "Geneid")

# tidy up sample names
names(all_201221) <- gsub(x = names(all_201221), pattern = "./", replacement = "")  
names(all_201221) <- gsub(x = names(all_201221), pattern = ".bam", replacement = "")
names(all_201221) <- gsub(x = names(all_201221), pattern = "1131-TP-1_CD8_NP16_", replacement = "")

# list of sample names
list_201221 <- as.data.frame(colnames(all_201221))
bulk_orf3a <- all_201203[,c(481:486,774:779, 1404:1408)]
bulk_np16_s34_m24 <- all_201221[,25:96]

list_bulk_np16_s34_m24 <- as.data.frame(colnames(bulk_np16_s34_m24))

bulk_data <- cbind(bulk_orf3a, bulk_np16_s34_m24)




############################

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
