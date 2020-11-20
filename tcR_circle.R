# combining mixcr parsed files into single dataframe 
# each row name is cell (some cells have multiple rows)

#---------------- CD4 TCR ---------------------------

library(tcR)

# TCR sequence database

data(genesegments)

# TCR
HUMAN_TRAV
HUMAN_TRAJ
HUMAN_TRBV
HUMAN_TRBD
HUMAN_TRBJ
HUMAN_TRBV_MITCR
HUMAN_TRBV_ALS
HUMAN_TRGV
HUMAN_TRGJ
HUMAN_TRDV
HUMAN_TRDD
HUMAN_TRDJ

# BCR
HUMAN_IGHV
HUMAN_IGHD
HUMAN_IGHJ
HUMAN_IGLV
HUMAN_IGLJ

# Nucleotide sequence and CDR3 position of each gene segment
genesegments$TRBV[1:10,]

library(data.table)


setwd('/t1-data/user/lfelce/TCR_analysis/cd4/')
listFiles = list.files()

# revise the structure of input files
for(i in 1:length(listFiles))
{
  DT = fread(listFiles[i])
  write.table(DT,listFiles[i],quote = F, row.names = F, sep = '\t')
}

#read files
mixcr <- parse.folder("/t1-data/user/lfelce/TCR_analysis/cd4/", 'mixcr')

all_metadata = fread('/t1-data/user/lfelce/TCR_analysis/metadata_cd4.txt', stringsAsFactors = F)
colnames(all_metadata)
all_metadata$patient <- as.character(all_metadata$patient)

# rename patients 022 and 025
i=1
for (i in (1:96)) {all_metadata[i,2] <- "022"}
for (i in (193:288)) {all_metadata[i,2] <- "025"}


# metadata has all samples - but not all samples have valid clones.

# names of files in folder- 398 files
name_list <- fread('/t1-data/user/lfelce/TCR_analysis/cd4_names.txt', stringsAsFactors = F, header=F)

# select only file names which have valid clones
metadata <- all_metadata[is.element(all_metadata$sample, name_list$V1),]names(mixcr)

mixcr_ordered = mixcr[order(names(mixcr))]
names(mixcr_ordered)

class(mixcr)

cd4_df <- do.call(rbind.data.frame, mixcr)

dim(cd4_df)

# new column for sample name

cd4_df$sample <- rownames(cd4_df)

# move sample name to 1st column

cd4_df <- cd4_df[,c(17, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16)]

# convert so can write data frame as .csv file

cd4_df <- apply(cd4_df,2,as.character)

# set working directory to export .csv

setwd('/t1-data/user/lfelce/TCR_analysis')

write.csv(cd4_df, file="all_cd4_tcr.csv", row.names=FALSE)



#----------------------- CD8 TCR ---------------------------

setwd('/t1-data/user/lfelce/TCR_analysis/cd8/')
listFiles = list.files()

# revise the structure of input files
for(i in 1:length(listFiles))
{
  DT = fread(listFiles[i])
  write.table(DT,listFiles[i],quote = F, row.names = F, sep = '\t')
}

#read files
mixcr <- parse.folder("/t1-data/user/lfelce/TCR_analysis/cd8/", 'mixcr')

all_metadata = fread('/t1-data/user/lfelce/TCR_analysis/metadata_cd8.txt', stringsAsFactors = F)
colnames(all_metadata)

# change 5 to 005
for (i in (1:192)) {all_metadata[i,2] <- "005"}

# metadata has all samples - but not all samples have valid clones.

# names of files in folder- 273 files
name_list <- fread('/t1-data/user/lfelce/TCR_analysis/cd8_names.txt', stringsAsFactors = F, header=F)

# select only file names which have valid clones
metadata <- all_metadata[is.element(all_metadata$sample, name_list$V1),]

cd8_df <- do.call(rbind.data.frame, mixcr)
dim(cd8_df)

# new column for sample name

cd8_df$sample <- rownames(cd8_df)

# move sample name to 1st column

cd8_df <- cd8_df[,c(17, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16)]

# convert so can write data frame as .csv file

cd8_df <- apply(cd8_df,2,as.character)

# set working directory to export .csv

setwd('/t1-data/user/lfelce/TCR_analysis')

write.csv(cd8_df, file="all_cd8_tcr.csv", row.names=FALSE)


#--------------------CD8-ORF samples------------------------

setwd('/t1-data/user/lfelce/TCR_analysis/cd8_orf/')
listFiles = list.files()

# revise the structure of input files
for(i in 1:length(listFiles))
{
  DT = fread(listFiles[i])
  write.table(DT,listFiles[i],quote = F, row.names = F, sep = '\t')
}

#read files
mixcr <- parse.folder("/t1-data/user/lfelce/TCR_analysis/cd8_orf/", 'mixcr')

all_metadata = fread('/t1-data/user/lfelce/TCR_analysis/metadata_cd8.txt', stringsAsFactors = F)
colnames(all_metadata)

# change 5 to 005
for (i in (1:192)) {all_metadata[i,2] <- "005"}

# metadata has all samples - but not all samples have valid clones.

# names of files in folder- 312 files
name_list <- fread('/t1-data/user/lfelce/TCR_analysis/cd8_orf_names.txt', stringsAsFactors = F, header=F)

# select only file names which have valid clones
metadata <- all_metadata[is.element(all_metadata$sample, name_list$V1),]

cd8_orf_df <- do.call(rbind.data.frame, mixcr)
dim(cd8_orf_df)

# new column for sample name

cd8_orf_df$sample <- rownames(cd8_orf_df)

# move sample name to 1st column

cd8_orf_df <- cd8_orf_df[,c(17, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16)]

# convert so can write data frame as .csv file

cd8_orf_df <- apply(cd8_orf_df,2,as.character)

# set working directory to export .csv

setwd('/t1-data/user/lfelce/TCR_analysis')

write.csv(cd8_orf_df, file="all_cd8_orf_tcr.csv", row.names=FALSE)

#---------------------Circular plot - CD4 ----------------------------------

library(circlize)
setwd('/t1-data/user/lfelce/TCR_analysis')

# separate S34 and M24 samples

cd4_s34 <- as.data.frame(rbind(cd4_df[c(1:134,282:400),]))

cd4_m24 <- as.data.frame(rbind(cd4_df[c(135:281,401:852),]))


# need to plot table with variable genes as rows and junction genes as columns

cd4_s34_table <- as.data.frame.matrix(table(cd4_s34$V.gene, cd4_s34$J.gene))
cd4_s34_table <- as.matrix(cd4_s34_table[-1,-1])

# alpha chain
cd4_s34_table_a <- as.data.frame.matrix(rbind(cd4_s34_table[c(1:31),]))
cd4_s34_table_a <- as.data.frame.matrix(cbind(cd4_s34_table_a[,c(1:36)]))
cd4_s34_table_a <- as.matrix(cd4_s34_table_a)

pdf('cd4_s34_a_chorddiagram.pdf', width = 12, height = 8, useDingbats = FALSE)
circos.clear()
set.seed(999)
chordDiagram(cd4_s34_table_a, annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.25, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)
dev.off()

# beta chain
cd4_s34_table_b <- as.data.frame.matrix(rbind(cd4_s34_table[c(32:72),]))
cd4_s34_table_b <- as.data.frame.matrix(cbind(cd4_s34_table_b[,c(37:47)]))
cd4_s34_table_b <- as.matrix(cd4_s34_table_b)

pdf('cd4_s34_b_chorddiagram.pdf', width = 12, height = 8, useDingbats = FALSE)
circos.clear()
set.seed(999)
chordDiagram(cd4_s34_table_b, annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.25, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)

dev.off()


# need to plot table with variable genes as rows and junction genes as columns

cd4_m24_table <- as.data.frame.matrix(table(cd4_m24$V.gene, cd4_m24$J.gene))
cd4_m24_table <- as.matrix(cd4_m24_table[-1,-1])

# alpha chain
cd4_m24_table_a <- as.data.frame.matrix(rbind(cd4_m24_table[c(1:46),]))
cd4_m24_table_a <- as.data.frame.matrix(cbind(cd4_m24_table_a[,c(1:49)]))
cd4_m24_table_a <- as.matrix(cd4_m24_table_a)

pdf('cd4_m24_a_chorddiagram.pdf', width = 12, height = 8, useDingbats = FALSE)
circos.clear()
set.seed(999)
chordDiagram(cd4_m24_table_a, annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.25, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)
dev.off()

# beta chain
cd4_m24_table_b <- as.data.frame.matrix(rbind(cd4_m24_table[c(47:110),]))
cd4_m24_table_b <- as.data.frame.matrix(cbind(cd4_m24_table_b[,c(50:63)]))
cd4_m24_table_b <- as.matrix(cd4_m24_table_b)

pdf('cd4_m24_b_chorddiagram.pdf', width = 12, height = 8, useDingbats = FALSE)
circos.clear()
set.seed(999)
chordDiagram(cd4_m24_table_b, annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.25, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)

dev.off()

#------------------- Circular plot CD8 tabulate and separate alpha and beta ------------------

cd8_df <- as.data.frame(cd8_df)
cd8_table <- as.data.frame.matrix(table(cd8_df$V.gene, cd8_df$J.gene))

# alpha chain
cd8_table_a <- as.data.frame.matrix(rbind(cd8_table[c(1:41),]))
cd8_table_a <- as.data.frame.matrix(cbind(cd8_table_a[,c(1:37)]))
cd8_table_a <- as.matrix(cd8_table_a)

pdf('cd8_a_chorddiagram.pdf', width = 12, height = 8, useDingbats = FALSE)
circos.clear()
set.seed(999)
chordDiagram(cd8_table_a, annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.25, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)
dev.off()

# beta chain
cd8_table_b <- as.data.frame.matrix(rbind(cd8_table[c(42:91),]))
cd8_table_b <- as.data.frame.matrix(cbind(cd8_table_b[,c(38:51)]))
cd8_table_b <- as.matrix(cd8_table_b)

pdf('cd8_b_chorddiagram.pdf', width = 12, height = 8, useDingbats = FALSE)
circos.clear()
set.seed(999)
chordDiagram(cd8_table_b, annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.25, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)

dev.off()

#----------------- Circular plot CD8 ORF --------------------------

cd8_orf_df <- as.data.frame(cd8_orf_df)
cd8_orf_table <- as.data.frame.matrix(table(cd8_orf_df$V.gene, cd8_orf_df$J.gene))
cd8_orf_table <- as.matrix(cd8_orf_table[-1,-1])

# alpha chain
cd8_orf_table_a <- as.data.frame.matrix(rbind(cd8_orf_table[c(1:48),]))
cd8_orf_table_a <- as.data.frame.matrix(cbind(cd8_orf_table_a[,c(1:48)]))
cd8_orf_table_a <- as.matrix(cd8_orf_table_a)

pdf('cd8_orf_a_chorddiagram.pdf', width = 12, height = 8, useDingbats = FALSE)
circos.clear()
set.seed(999)
chordDiagram(cd8_orf_table_a, annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.25, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)
dev.off()

# beta chain
cd8_orf_table_b <- as.data.frame.matrix(rbind(cd8_orf_table[c(49:160),]))
cd8_orf_table_b <- as.data.frame.matrix(cbind(cd8_orf_table_b[,c(49:64)]))
cd8_orf_table_b <- as.matrix(cd8_orf_table_b)

pdf('cd8_orf_b_chorddiagram.pdf', width = 12, height = 8, useDingbats = FALSE)
circos.clear()
set.seed(999)
chordDiagram(cd8_orf_table_b, annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.25, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)

dev.off()

#-----------------------Separate by patient ------------------------

# CD8 005

cd8_005 <- cd8_df[1:170,]
cd8_005_table <- as.data.frame.matrix(table(cd8_005$V.gene, cd8_005$J.gene))
cd8_005_table <- as.matrix(cd8_005_table)
circos.clear()
set.seed(999)
chordDiagram(cd8_005_table, annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.25, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)


# CD8 1131TP1
cd8_1131TP1 <- cd8_df[171:175,]
cd8_1131TP1_table <- as.data.frame.matrix(table(cd8_1131TP1$V.gene, cd8_1131TP1$J.gene))
cd8_1131TP1_table <- as.matrix(cd8_1131TP1_table)
circos.clear()
set.seed(999)
chordDiagram(cd8_1131TP1_table, annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.25, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)

# CD8 1131TP2
cd8_1131TP2 <- cd8_df[176:298,]
cd8_1131TP2_table <- as.data.frame.matrix(table(cd8_1131TP2$V.gene, cd8_1131TP2$J.gene))
cd8_1131TP2_table <- as.matrix(cd8_1131TP2_table)
circos.clear()
set.seed(999)
chordDiagram(cd8_1131TP2_table, annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.25, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)

# CD8 1153
cd8_1153 <- cd8_df[299:407,]
cd8_1153_table <- as.data.frame.matrix(table(cd8_1153$V.gene, cd8_1153$J.gene))
cd8_1153_table <- as.matrix(cd8_1153_table)
circos.clear()
set.seed(999)
chordDiagram(cd8_1153_table, annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.25, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)

# CD8 1201
cd8_1201 <- cd8_df[408:563,]
cd8_1201_table <- as.data.frame.matrix(table(cd8_1201$V.gene, cd8_1201$J.gene))
cd8_1201_table <- as.matrix(cd8_1201_table)
circos.clear()
set.seed(999)
chordDiagram(cd8_1201_table, annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.25, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)

# CD4 022
cd4_022 <- cd4_df[1:134,]
cd4_022_table <- as.data.frame.matrix(table(cd4_022$V.gene, cd4_022$J.gene))
cd4_022_table <- as.matrix(cd4_022_table)
circos.clear()
set.seed(999)
chordDiagram(cd4_022_table, annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.25, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)

# CD4 025
cd4_025 <- cd4_df[135:281,]
cd4_025_table <- as.data.frame.matrix(table(cd4_025$V.gene, cd4_025$J.gene))
cd4_025_table <- as.matrix(cd4_025_table)
circos.clear()
set.seed(999)
chordDiagram(cd4_025_table, annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.25, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)

# CD4 1062
cd4_1062 <- cd4_df[282:400,]
cd4_1062_table <- as.data.frame.matrix(table(cd4_1062$V.gene, cd4_1062$J.gene))
cd4_1062_table <- as.matrix(cd4_1062_table)
circos.clear()
set.seed(999)
chordDiagram(cd4_1062_table, annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.25, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)

# CD4 1493
cd4_1493 <- cd4_df[401:555,]
cd4_1493_table <- as.data.frame.matrix(table(cd4_1493$V.gene, cd4_1493$J.gene))
cd4_1493_table <- as.matrix(cd4_1493_table)
# if not clear from chord diagram, then view table and search for numbers like 2, 3, 4, etc to see highest hits
circos.clear()
set.seed(999)
chordDiagram(cd4_1493_table, annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.25, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)

# CD4 1504
cd4_1504 <- cd4_df[556:709,]
cd4_1504_table <- as.data.frame.matrix(table(cd4_1504$V.gene, cd4_1504$J.gene))
cd4_1504_table <- as.matrix(cd4_1504_table)
# if not clear from chord diagram, then view table and search for numbers like 2, 3, 4, etc to see highest hits
# copy and paste into Excel as can freeze 1st column to see V genes better
# and can do conditional formatting
circos.clear()
set.seed(999)
chordDiagram(cd4_1504_table, annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.25, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)

# CD4 1525
cd4_1525 <- cd4_df[710:852,]
cd4_1525_table <- as.data.frame.matrix(table(cd4_1525$V.gene, cd4_1525$J.gene))
cd4_1525_table <- as.matrix(cd4_1525_table)
# if not clear from chord diagram, then view table and search for numbers like 2, 3, 4, etc to see highest hits
# copy and paste into Excel as can freeze 1st column to see V genes better
# and can do conditional formatting
circos.clear()
set.seed(999)
chordDiagram(cd4_1525_table, annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.25, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)







#------------------------- Finding alpha and beta dominant pair for each patient --------------------
library(dplyr)

## C8 NP-16 patients
# 005
datalist = list()

for (i in (1:91)) {
  # ... make some data
  dat <- data.frame(c(mixcr[[i]][7], mixcr[[i]][8]))
  dat$i <- i # maybe you want to keep track of which iteration produced it?
  datalist[[i]] <- dat # add it to your list
}

big_data = do.call(rbind, datalist)
big_data1 <- big_data %>% group_by(i) %>% filter(n() >= 2)

big_data2 <- mutate(big_data1, VJ = paste(V.gene, J.gene))
# copy and paste into Excel, sort manually into Alpha and Beta columns
# copy into Notepad, copy into cluster Excel and save as .csv

pairs_005 <- read.csv("pairs_005.csv")
pairs_005_table <- as.data.frame.matrix(table(pairs_005$Beta, pairs_005$Alpha))
pairs_005_table <- as.matrix(pairs_005_table)

circos.clear()
set.seed(999)
chordDiagram(pairs_005_table, annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.25, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)

# 1131-TP-1 - only 1 suitable pairing, so can just check table by eye!
datalist = list()

for (i in (92:95)) {
  # ... make some data
  dat <- data.frame(c(mixcr[[i]][7], mixcr[[i]][8]))
  dat$i <- i # maybe you want to keep track of which iteration produced it?
  datalist[[i]] <- dat # add it to your list
}

big_data = do.call(rbind, datalist)

# 1131-TP-2
datalist = list()

for (i in (96:155)) {
  # ... make some data
  dat <- data.frame(c(mixcr[[i]][7], mixcr[[i]][8]))
  dat$i <- i # maybe you want to keep track of which iteration produced it?
  datalist[[i]] <- dat # add it to your list
}

big_data = do.call(rbind, datalist)
big_data1 <- big_data %>% group_by(i) %>% filter(n() >= 2)

big_data2 <- mutate(big_data1, VJ = paste(V.gene, J.gene))

pairs <- read.csv("pairs_1131-TP-2.csv")
pairs_table <- as.data.frame.matrix(table(pairs$Beta, pairs$Alpha))
pairs_table <- as.matrix(pairs_table)

circos.clear()
set.seed(999)
chordDiagram(pairs_table, annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.25, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)

# 1153
datalist = list()

for (i in (156:206)) {
  # ... make some data
  dat <- data.frame(c(mixcr[[i]][7], mixcr[[i]][8]))
  dat$i <- i # maybe you want to keep track of which iteration produced it?
  datalist[[i]] <- dat # add it to your list
}

big_data = do.call(rbind, datalist)
big_data1 <- big_data %>% group_by(i) %>% filter(n() >= 2)

big_data2 <- mutate(big_data1, VJ = paste(V.gene, J.gene))

pairs <- read.csv("pairs_1153.csv")
pairs_table <- as.data.frame.matrix(table(pairs$Beta, pairs$Alpha))
pairs_table <- as.matrix(pairs_table)

circos.clear()
set.seed(999)
chordDiagram(pairs_table, annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.25, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)

# 1201-TP-2
datalist = list()

for (i in (207:273)) {
  # ... make some data
  dat <- data.frame(c(mixcr[[i]][7], mixcr[[i]][8]))
  dat$i <- i # maybe you want to keep track of which iteration produced it?
  datalist[[i]] <- dat # add it to your list
}

big_data = do.call(rbind, datalist)
big_data1 <- big_data %>% group_by(i) %>% filter(n() >= 2)

big_data2 <- mutate(big_data1, VJ = paste(V.gene, J.gene))

pairs <- read.csv("pairs_1201-TP-2.csv")
pairs_table <- as.data.frame.matrix(table(pairs$Beta, pairs$Alpha))
pairs_table <- as.matrix(pairs_table)

circos.clear()
set.seed(999)
chordDiagram(pairs_table, annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.25, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)
