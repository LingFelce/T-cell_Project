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



#----------------------- CD8 TCR

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


#---------------------Circular plot - CD4

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

# cd8 tabulate and separate alpha and beta

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


