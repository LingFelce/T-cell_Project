##### Updated MiXCR-SmartSeq.sh output####
# Isar modified MiXCR script so that .txt output for each cell is TRA.txt, TRB.txt, TRD.txt, TRG.txt file
library(tcR)

#-------------  List of human of TCR and Ig gene segments

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

#------------- CD8 NP16 -------------------------

setwd('/t1-data/user/lfelce/TCR_analysis/cd8_np16_new/')
listFiles = list.files()

# revise the structure of input files
for(i in 1:length(listFiles))
{
  DT = fread(listFiles[i])
  write.table(DT,listFiles[i],quote = F, row.names = F, sep = '\t')
}

library(stringr)

# read in TRA names and convert to list
tra_names <- fread('/t1-data/user/lfelce/TCR_analysis/cd8_np16_tra_names.txt', stringsAsFactors = F, header=F)

tra_names <- as.list(as.data.frame(t(tra_names)))

# remove ./ at start of name and replace with file path
tra_names <- tra_names %>% str_replace("./*", "")

tra_names <- paste("/t1-data/user/lfelce/TCR_analysis/cd8_np16_new/", tra_names, sep="")

# read in TRB names and convert to list
trb_names <- fread('/t1-data/user/lfelce/TCR_analysis/cd8_np16_trb_names.txt', stringsAsFactors = F, header=F)

trb_names <- as.list(as.data.frame(t(trb_names)))

# remove ./ at start of name and replace with file path
trb_names <- trb_names %>% str_replace("./*", "")

trb_names <- paste("/t1-data/user/lfelce/TCR_analysis/cd8_np16_new/", trb_names, sep="")

# parse mixcr files
mixcr_a <- parse.file.list(tra_names, "mixcr")
mixcr_b <- parse.file.list(trb_names, "mixcr")

# sort alphabetically
mixcr_a <- mixcr_a[order(names(mixcr_a))]
mixcr_b <- mixcr_b[order(names(mixcr_b))]

# create list of cell numbers and cell names
# different lengths so same cell will be different number in a or b
mixcr_a_names <- as.data.frame(names(mixcr_a))
mixcr_b_names <- as.data.frame(names(mixcr_b))

mixcr_a_names <-tibble::rownames_to_column(mixcr_a_names, "cell_number")
mixcr_b_names <- tibble::rownames_to_column(mixcr_b_names, "cell_number")

# rename columns
colnames(mixcr_a_names) <- c("cell_number", "cell_name")
colnames(mixcr_b_names) <- c("cell_number", "cell_name")

# convert mixcr lists to dataframe with just V.gene and J.gene info

# mixcr_a
datalist = list()
for (i in (1:length(mixcr_a))) {
  dat <- data.frame(c(mixcr_a[[i]][7], mixcr_a[[i]][8]))
  dat$i <- i # keep track of which iteration produced it
  datalist[[i]] <- dat # add it to list
}
# combine columns for each cell, select only cells with only 2 rows (dual alpha)
big_data = do.call(rbind, datalist)
tra <- big_data %>% group_by(i) %>% filter(n() <= 2)
colnames(tra) <- c("TRAV", "TRAJ", "cell_number")

# mixcr_b
datalist = list()
for (i in (1:length(mixcr_b))) {
  dat <- data.frame(c(mixcr_b[[i]][7], mixcr_b[[i]][8]))
  dat$i <- i # keep track of which iteration produced it
  datalist[[i]] <- dat # add it to list
}
# combine columns for each cell, select only cells with only 1 row (single beta)
big_data = do.call(rbind, datalist)
trb <- big_data %>% group_by(i) %>% filter(n() == 1)
colnames(trb) <- c("TRBV", "TRBJ", "cell_number")

# combine TRA and TRB dataframes
cd8_np16 <- merge(tra,trb, by="cell_number")

# create list of cell numbers and cell names
mixcr_a_names <-tibble::rownames_to_column(mixcr_a_names, "cell_number")
mixcr_b_names <- tibble::rownames_to_column(mixcr_b_names, "cell_number")

# merge names with main dataframe
cd8_np16 <- merge(cd8_np16, mixcr_b_names, by="cell_number")

cd8_np16 <- mutate(cd8_np16, alpha=paste(TRAV, TRAJ, sep="_"))

cd8_np16 <- mutate(cd8_np16, beta=paste(TRBV, TRBJ, sep="_"))

# tabulate to get dominant alpha-beta pairing
# 005 
# 1131-TP-1 - no results
# 1131-TP-2 
# 1153 
# 1201-TP-2 
library(circlize)

patient_005 <- cd8_np16[cd8_np16$`names(mixcr_b)` %like% "005",]
cd8_np16_005 <- as.matrix(as.data.frame.matrix(table(patient_005$alpha, patient_005$beta)))

patient_1131 <- cd8_np16[cd8_np16$`names(mixcr_b)` %like% "1131-TP-2",]
cd8_np16_1131 <- as.matrix(as.data.frame.matrix(table(patient_1131$alpha, patient_1131$beta)))

patient_1153 <- cd8_np16[cd8_np16$`names(mixcr_b)` %like% "1153",]
cd8_np16_1153 <- as.matrix(as.data.frame.matrix(table(patient_1153$alpha, patient_1153$beta)))

patient_1201 <- cd8_np16[cd8_np16$`names(mixcr_b)` %like% "1201-TP-2",]
cd8_np16_1201 <- as.matrix(as.data.frame.matrix(table(patient_1201$alpha, patient_1201$beta)))

setwd("/t1-data/user/lfelce/TCR_analysis/new_mixcr_results/")

list <- c("cd8_np16_005","cd8_np16_1131", "cd8_np16_1153", "cd8_np16_1201")

for (i in 1:length(list)) {
  circos.clear()
  set.seed(999)
  pdf(paste((list[i]), "_chorddiagram.pdf", sep=""), width = 16, height = 12, useDingbats = FALSE)
  chordDiagram(get(list[i]), annotationTrack = "grid", preAllocateTracks = 1)
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
    circos.axis(h = "top", labels.cex = 0.25, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
  }, bg.border = NA)
  dev.off()
}

#------------- CD8 ORF3a-28 -------------------------

setwd('/t1-data/user/lfelce/TCR_analysis/cd8_orf_new/')
listFiles = list.files()

# revise the structure of input files
for(i in 1:length(listFiles))
{
  DT = fread(listFiles[i])
  write.table(DT,listFiles[i],quote = F, row.names = F, sep = '\t')
}

library(stringr)

# read in TRA names and convert to list
tra_names <- fread('/t1-data/user/lfelce/TCR_analysis/cd8_orf_tra_names.txt', stringsAsFactors = F, header=F)

tra_names <- as.list(as.data.frame(t(tra_names)))

# remove ./ at start of name and replace with file path
tra_names <- tra_names %>% str_replace("./*", "")

tra_names <- paste("/t1-data/user/lfelce/TCR_analysis/cd8_orf_new/", tra_names, sep="")

# read in TRB names and convert to list
trb_names <- fread('/t1-data/user/lfelce/TCR_analysis/cd8_orf_trb_names.txt', stringsAsFactors = F, header=F)

trb_names <- as.list(as.data.frame(t(trb_names)))

# remove ./ at start of name and replace with file path
trb_names <- trb_names %>% str_replace("./*", "")

trb_names <- paste("/t1-data/user/lfelce/TCR_analysis/cd8_orf_new/", trb_names, sep="")

# parse mixcr files
mixcr_a <- parse.file.list(tra_names, "mixcr")
mixcr_b <- parse.file.list(trb_names, "mixcr")

# sort alphabetically
mixcr_a <- mixcr_a[order(names(mixcr_a))]
mixcr_b <- mixcr_b[order(names(mixcr_b))]

# create list of cell numbers and cell names
# different lengths so same cell will be different number in a or b
mixcr_a_names <- as.data.frame(names(mixcr_a))
mixcr_b_names <- as.data.frame(names(mixcr_b))

mixcr_a_names <-tibble::rownames_to_column(mixcr_a_names, "cell_number")
mixcr_b_names <- tibble::rownames_to_column(mixcr_b_names, "cell_number")

# rename columns
colnames(mixcr_a_names) <- c("cell_number", "cell_name")
colnames(mixcr_b_names) <- c("cell_number", "cell_name")

# convert mixcr lists to dataframe with just V.gene and J.gene info

# mixcr_a
datalist = list()
for (i in (1:(length(mixcr_a)))) {
  dat <- data.frame(c(mixcr_a[[i]][7], mixcr_a[[i]][8]))
  dat$i <- i # keep track of which iteration produced it
  datalist[[i]] <- dat # add it to list
}
# combine columns for each cell, select only cells with only 2 rows (dual receptor)
big_data = do.call(rbind, datalist)
tra <- big_data %>% group_by(i) %>% filter(n() <= 2)
colnames(tra) <- c("TRAV", "TRAJ", "cell_number")
tra <- merge(tra, mixcr_a_names, by="cell_number")


# mixcr_b
datalist = list()
for (i in (1:(length(mixcr_b)))) {
  dat <- data.frame(c(mixcr_b[[i]][7], mixcr_b[[i]][8]))
  dat$i <- i # keep track of which iteration produced it
  datalist[[i]] <- dat # add it to list
}
# combine columns for each cell, select only cells with only 1 row (single receptor)
big_data = do.call(rbind, datalist)
trb <- big_data %>% group_by(i) %>% filter(n() == 1)
colnames(trb) <- c("TRBV", "TRBJ", "cell_number")
trb <- merge(trb, mixcr_b_names, by="cell_number")

# combine TRA and TRB dataframes
cd8_orf <- merge(tra,trb, by="cell_name")


# merge names with main dataframe
cd8_orf <- merge(cd8_orf, mixcr_b_names, by="cell_number")

cd8_orf <- mutate(cd8_orf, alpha=paste(TRAV, TRAJ, sep="_"))

cd8_orf <- mutate(cd8_orf, beta=paste(TRBV, TRBJ, sep="_"))

# tabulate to get dominant alpha-beta pairing
# 1105 rows 
# 1134-TP-2 
# 1525-TP-1 
# 1525-TP-2

library(circlize)

patient_1105 <- cd8_orf[cd8_orf$cell_name %like% "1105", ]
cd8_orf_1105 <- as.matrix(as.data.frame.matrix(table(patient_1105$alpha, patient_1105$beta)))

patient_1134 <- cd8_orf[cd8_orf$cell_name %like% "1134-TP-2", ]
cd8_orf_1134 <- as.matrix(as.data.frame.matrix(table(patient_1134$alpha, patient_1134$beta)))

patient_1525TP1 <- cd8_orf[cd8_orf$cell_name %like% "1525-TP-1", ]
cd8_orf_1525TP1 <- as.matrix(as.data.frame.matrix(table(patient_1525TP1$alpha, patient_1525TP1$beta)))

patient_1525TP2 <- cd8_orf[cd8_orf$cell_name %like% "1525-TP-2", ]
cd8_orf_1525TP2 <- as.matrix(as.data.frame.matrix(table(patient_1525TP2$alpha, patient_1525TP2$beta)))

setwd("/t1-data/user/lfelce/TCR_analysis/new_mixcr_results/")

list <- c("cd8_orf_1105","cd8_orf_1134", "cd8_orf_1525TP1", "cd8_orf_1525TP2")

for (i in 1:length(list)) {
  circos.clear()
  set.seed(999)
  pdf(paste((list[i]), "_chorddiagram.pdf", sep=""), width = 16, height = 12, useDingbats = FALSE)
  chordDiagram(get(list[i]), annotationTrack = "grid", preAllocateTracks = 1)
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
    circos.axis(h = "top", labels.cex = 0.25, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
  }, bg.border = NA)
  dev.off()
}

#------------- CD4 S34 and M24 -------------------------

setwd('/t1-data/user/lfelce/TCR_analysis/cd4_new/')
listFiles = list.files()

# revise the structure of input files
for(i in 1:length(listFiles))
{
  DT = fread(listFiles[i])
  write.table(DT,listFiles[i],quote = F, row.names = F, sep = '\t')
}

library(stringr)

# read in TRA names and convert to list
tra_names <- fread('/t1-data/user/lfelce/TCR_analysis/cd4_tra_names.txt', stringsAsFactors = F, header=F)

tra_names <- as.list(as.data.frame(t(tra_names)))

# remove ./ at start of name and replace with file path
tra_names <- tra_names %>% str_replace("./*", "")

tra_names <- paste("/t1-data/user/lfelce/TCR_analysis/cd4_new/", tra_names, sep="")

# read in TRB names and convert to list
trb_names <- fread('/t1-data/user/lfelce/TCR_analysis/cd4_trb_names.txt', stringsAsFactors = F, header=F)

trb_names <- as.list(as.data.frame(t(trb_names)))

# remove ./ at start of name and replace with file path
trb_names <- trb_names %>% str_replace("./*", "")

trb_names <- paste("/t1-data/user/lfelce/TCR_analysis/cd4_new/", trb_names, sep="")

# parse mixcr files
mixcr_a <- parse.file.list(tra_names, "mixcr")
mixcr_b <- parse.file.list(trb_names, "mixcr")

# sort alphabetically
mixcr_a <- mixcr_a[order(names(mixcr_a))]
mixcr_b <- mixcr_b[order(names(mixcr_b))]

# convert mixcr lists to dataframe with just V.gene and J.gene info

# mixcr_a
mixcr_a_names <- as.data.frame(names(mixcr_a))

datalist = list()
for (i in (1:(length(mixcr_a)))) {
  dat <- data.frame(c(mixcr_a[[i]][7], mixcr_a[[i]][8]))
  dat$i <- i # keep track of which iteration produced it
  datalist[[i]] <- dat # add it to list
}
# combine columns for each cell, select only cells with only 2 rows (dual alpha)
big_data = do.call(rbind, datalist)
tra <- big_data %>% group_by(i) %>% filter(n() <= 2)
colnames(tra) <- c("TRAV", "TRAJ", "cell_number")

# mixcr_b
mixcr_b_names <- as.data.frame(names(mixcr_b))

datalist = list()
for (i in (1:(length(mixcr_b)))) {
  dat <- data.frame(c(mixcr_b[[i]][7], mixcr_b[[i]][8]))
  dat$i <- i # keep track of which iteration produced it
  datalist[[i]] <- dat # add it to list
}
# combine columns for each cell, select only cells with only 1 row (single beta)
big_data = do.call(rbind, datalist)
trb <- big_data %>% group_by(i) %>% filter(n() == 1)
colnames(trb) <- c("TRBV", "TRBJ", "cell_number")

# combine TRA and TRB dataframes
cd4 <- merge(tra,trb, by="cell_number")

# create list of cell numbers and cell names
mixcr_a_names <-tibble::rownames_to_column(mixcr_a_names, "cell_number")
mixcr_b_names <- tibble::rownames_to_column(mixcr_b_names, "cell_number")

# merge names with main dataframe
cd4 <- merge(cd4, mixcr_b_names, by="cell_number")

cd4 <- mutate(cd4, alpha=paste(TRAV, TRAJ, sep="_"))

cd4 <- mutate(cd4, beta=paste(TRBV, TRBJ, sep="_"))

# tabulate to get dominant alpha-beta pairing
# 022 rows 
# 025
# 1062
# 1493
# 1504
# 1525

library(circlize)

patient_022 <- cd4[cd4$`names(mixcr_b)` %like% "022", ]
cd4_s34_022 <- as.matrix(as.data.frame.matrix(table(patient_022$alpha, patient_022$beta)))

patient_025 <- cd4[cd4$`names(mixcr_b)` %like% "025", ]
cd4_m24_025 <- as.matrix(as.data.frame.matrix(table(patient_025$alpha, patient_025$beta)))

patient_1062 <- cd4[cd4$`names(mixcr_b)` %like% "1062", ]
cd4_s34_1062 <- as.matrix(as.data.frame.matrix(table(patient_1062$alpha, patient_1062$beta)))

patient_1493 <- cd4[cd4$`names(mixcr_b)` %like% "1493", ]
cd4_m24_1493 <- as.matrix(as.data.frame.matrix(table(patient_1493$alpha, patient_1493$beta)))

patient_1504 <- cd4[cd4$`names(mixcr_b)` %like% "1504", ]
cd4_m24_1504 <- as.matrix(as.data.frame.matrix(table(patient_1504$alpha, patient_1504$beta)))

patient_1525 <- cd4[cd4$`names(mixcr_b)` %like% "1525", ]
cd4_m24_1525 <- as.matrix(as.data.frame.matrix(table(patient_1525$alpha, patient_1525$beta)))

setwd("/t1-data/user/lfelce/TCR_analysis/new_mixcr_results/")

list <- c("cd4_s34_022", "cd4_m24_025", "cd4_s34_1062", "cd4_m24_1493", "cd4_m24_1504", "cd4_m24_1525")

for (i in 1:length(list)) {
  circos.clear()
  set.seed(999)
  pdf(paste((list[i]), "_chorddiagram.pdf", sep=""), width = 16, height = 12, useDingbats = FALSE)
  chordDiagram(get(list[i]), annotationTrack = "grid", preAllocateTracks = 1)
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
    circos.axis(h = "top", labels.cex = 0.25, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
  }, bg.border = NA)
  dev.off()
}
