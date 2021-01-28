##### Updated MiXCR-SmartSeq.sh output####
# Isar modified MiXCR script so that .txt output for each cell is TRA.txt, TRB.txt, TRD.txt, TRG.txt file
library(tcR)

#-------------  List of human of TCR and Ig gene segments ----------------------

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

# create list of file names
tra_list <- list.files(path = ".", recursive = TRUE,
                       pattern = "\\TRA.txt$", 
                       full.names = TRUE)
tra_list <- tra_list %>% str_replace("./*", "")

trb_list <- list.files(path = ".", recursive = TRUE,
                       pattern = "\\TRB.txt$", 
                       full.names = TRUE)
trb_list <- trb_list %>% str_replace("./*", "")


# parse in files with tcR
mixcr_a <- parse.file.list(tra_list, "mixcr")
mixcr_b <- parse.file.list(trb_list, "mixcr")

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

# convert mixcr lists to dataframe with Read count and proportion and V and J gene info
# re-order rows for each cell based on highest clone count first
# take top 2 rows to merge into large dataframe with all cells

# mixcr_a
datalist = list()
for (i in (1:length(mixcr_a))) {
  dat <- data.frame(c(mixcr_a[[i]][3], mixcr_a[[i]][4],mixcr_a[[i]][7], mixcr_a[[i]][8]))
  dat$i <- i # keep track of which iteration produced it
  dat <- dat[order(dat$Read.count),]
  dat <- dat[1:2,]
  datalist[[i]] <- dat # add it to list
}
# combine columns for each cell and filter any low false 2nd alphas
big_data = do.call(rbind, datalist)
tra <- big_data 
colnames(tra) <- c("clone_count", "clone_fraction", "TRAV", "TRAJ", "cell_number")
tra <- merge(tra, mixcr_a_names, by="cell_number")
tra <- tra[tra$clone_fraction > 0.01,]


# mixcr_b
datalist = list()
for (i in (1:length(mixcr_b))) {
  dat <- data.frame(c(mixcr_b[[i]][3], mixcr_b[[i]][4], mixcr_b[[i]][7], mixcr_b[[i]][8]))
  dat$i <- i # keep track of which iteration produced it
  dat <- dat[order(-dat$Read.count),]
  dat <- dat[1,]
  datalist[[i]] <- dat # add it to list
}
# combine columns for each cell
big_data = do.call(rbind, datalist)
trb <- big_data 
colnames(trb) <- c("clone_count", "clone_fraction", "TRBV", "TRBJ", "cell_number")
trb <- merge(trb, mixcr_b_names, by="cell_number")


# combine TRA and TRB dataframes - some cells will have only 1 alpha or 1 beta, so keep all rows
cd8_np16 <- merge(tra,trb, by="cell_name", all= TRUE)

cd8_np16 <- mutate(cd8_np16, alpha=paste(TRAV, TRAJ, sep="_"))

cd8_np16 <- mutate(cd8_np16, beta=paste(TRBV, TRBJ, sep="_"))

setwd('/t1-data/user/lfelce/TCR_analysis/new_mixcr_results/')

write.csv(cd8_np16, "cd8_np16_all_sc_tcr.csv")

# tabulate to get dominant alpha-beta pairing
# 005 
# 1131-TP-1 - no results
# 1131-TP-2 
# 1153 
# 1201-TP-2 
library(circlize)

patient_005 <- cd8_np16[cd8_np16$cell_name %like% "005",]
cd8_np16_005 <- as.matrix(as.data.frame.matrix(table(patient_005$alpha, patient_005$beta)))

patient_1131 <- cd8_np16[cd8_np16$cell_name %like% "1131-TP-2",]
cd8_np16_1131 <- as.matrix(as.data.frame.matrix(table(patient_1131$alpha, patient_1131$beta)))

patient_1153 <- cd8_np16[cd8_np16$cell_name %like% "1153",]
cd8_np16_1153 <- as.matrix(as.data.frame.matrix(table(patient_1153$alpha, patient_1153$beta)))

patient_1201 <- cd8_np16[cd8_np16$cell_name %like% "1201-TP-2",]
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

cd8_orf <- mutate(cd8_orf, alpha=paste(TRAV, TRAJ, sep="_"))

cd8_orf <- mutate(cd8_orf, beta=paste(TRBV, TRBJ, sep="_"))

setwd('/t1-data/user/lfelce/TCR_analysis/new_mixcr_results/')

write.csv(cd8_orf, "cd8_orf_sc_tcr.csv")

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
# combine columns for each cell, select only cells with only 2 rows (dual alpha)
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
# combine columns for each cell, select only cells with only 1 row (single beta)
big_data = do.call(rbind, datalist)
trb <- big_data %>% group_by(i) %>% filter(n() == 1)
colnames(trb) <- c("TRBV", "TRBJ", "cell_number")
trb <- merge(trb, mixcr_b_names, by="cell_number")

# combine TRA and TRB dataframes  
cd4 <- merge(tra,trb, by="cell_name")

cd4 <- mutate(cd4, alpha=paste(TRAV, TRAJ, sep="_"))

cd4 <- mutate(cd4, beta=paste(TRBV, TRBJ, sep="_"))

setwd('/t1-data/user/lfelce/TCR_analysis/new_mixcr_results/')

write.csv(cd4, "cd4_sc_tcr.csv")


# tabulate to get dominant alpha-beta pairing
library(circlize)

# combined cd4 patients

cd4_freq <- as.data.frame(table(cd4$alpha, cd4$beta))
cd4_all <- as.matrix(as.data.frame.matrix(table(cd4$alpha, cd4$beta)))

setwd('/t1-data/user/lfelce/TCR_analysis/new_mixcr_results')

circos.clear()
set.seed(999)
pdf('cd4_all_chorddiagram.pdf', width = 18, height = 14, useDingbats = FALSE)
chordDiagram(cd4_all, annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.25, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)
dev.off()

# combine only non-unique pairings for all cd4

combined_cd4 <- merge(cd4_s34_freq, cd4_m24_freq, by=c("Var1", "Var2"), all=T)
combined_cd4 <- na.omit(combined_cd4)

# look at alpha and beta chains separately - alpha chain first
tra_s34 <- tra[c(1:53, 113:159),]
tra_m24 <- tra[-c(1:53, 113:159),]

tra_s34_freq <- as.data.frame(table(tra_s34$TRAV, tra_s34$TRAJ))
tra_m24_freq <- as.data.frame(table(tra_m24$TRAV, tra_m24$TRAJ))

combined_tra <- merge(tra_s34_freq, tra_m24_freq, by=c("Var1", "Var2"), all = T)
combined_tra[combined_tra == 0] <- NA
combined_tra <- na.omit(combined_tra)

# look at alpha and beta chains separately - beta chain next
trb_s34 <- trb[c(1:29, 73:93),]
trb_m24 <- trb[-c(1:29, 73:93),]

trb_s34_freq <- as.data.frame(table(trb_s34$TRBV, trb_s34$TRBJ))
trb_m24_freq <- as.data.frame(table(trb_m24$TRBV, trb_m24$TRBJ))

combined_trb <- merge(trb_s34_freq, trb_m24_freq, by=c("Var1", "Var2"), all = T)
combined_trb[combined_trb == 0] <- NA
combined_trb <- na.omit(combined_trb)

# separate circle plots for each patient
# 022  
# 025
# 1062
# 1493
# 1504
# 1525

patient_022 <- cd4[cd4$cell_name %like% "022", ]
cd4_s34_022 <- as.matrix(as.data.frame.matrix(table(patient_022$alpha, patient_022$beta)))

patient_025 <- cd4[cd4$cell_name %like% "025", ]
cd4_m24_025 <- as.matrix(as.data.frame.matrix(table(patient_025$alpha, patient_025$beta)))

patient_1062 <- cd4[cd4$cell_name %like% "1062", ]
cd4_s34_1062 <- as.matrix(as.data.frame.matrix(table(patient_1062$alpha, patient_1062$beta)))

patient_1493 <- cd4[cd4$cell_name %like% "1493", ]
cd4_m24_1493 <- as.matrix(as.data.frame.matrix(table(patient_1493$alpha, patient_1493$beta)))

patient_1504 <- cd4[cd4$cell_name %like% "1504", ]
cd4_m24_1504 <- as.matrix(as.data.frame.matrix(table(patient_1504$alpha, patient_1504$beta)))

patient_1525 <- cd4[cd4$cell_name %like% "1525", ]
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


#-------------------------- CD4 S34 patients -----------------------

cd4_s34 <- cd4[c(1:23, 61:79),]

cd4_s34_freq <- as.data.frame(table(cd4_s34$alpha, cd4_s34$beta))
cd4_s34_all <- as.matrix(as.data.frame.matrix(table(cd4_s34$alpha, cd4_s34$beta)))

setwd('/t1-data/user/lfelce/TCR_analysis/new_mixcr_results')

circos.clear()
set.seed(999)
pdf('cd4_s34_all_chorddiagram.pdf', width = 18, height = 14, useDingbats = FALSE)
chordDiagram(cd4_s34_all, annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.25, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)
dev.off()


#--------------------------- CD4 M24 patients ----------------------

# individually CD4 M24 patients don't appear to have dominant alpha-beta pair
# look at CD4 M24 patients altogether and see if when combined there appears to be a dominant pair

cd4_m24 <- cd4[-c(1:23, 61:79),]
cd4_m24_freq <- as.data.frame(table(cd4_m24$alpha, cd4_m24$beta))
cd4_m24_all <- as.matrix(as.data.frame.matrix(table(cd4_m24$alpha, cd4_m24$beta)))

setwd('/t1-data/user/lfelce/TCR_analysis/new_mixcr_results')

circos.clear()
set.seed(999)
pdf('cd4_m24_all_chorddiagram.pdf', width = 16, height = 12, useDingbats = FALSE)
chordDiagram(cd4_m24_all, annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.25, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)
dev.off()

# separate by mild and severe

cd4_m24_mild <- cd4_m24[1:61,]
cd4_m24_mild_freq <- as.data.frame(table(cd4_m24_mild$alpha, cd4_m24_mild$beta))
cd4_m24_mild_table <- as.matrix(as.data.frame.matrix(table(cd4_m24_mild$alpha, cd4_m24_mild$beta)))

circos.clear()
set.seed(999)
pdf('cd4_m24_mild_chorddiagram.pdf', width = 16, height = 12, useDingbats = FALSE)
chordDiagram(cd4_m24_mild_table, annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.25, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)
dev.off()


cd4_m24_severe <- cd4_m24[62:95,]
cd4_m24_severe_freq <- as.data.frame(table(cd4_m24_severe$alpha, cd4_m24_severe$beta))
cd4_m24_severe_table <- as.matrix(as.data.frame.matrix(table(cd4_m24_severe$alpha, cd4_m24_severe$beta)))

circos.clear()
set.seed(999)
pdf('cd4_m24_severe_chorddiagram.pdf', width = 16, height = 12, useDingbats = FALSE)
chordDiagram(cd4_m24_severe_table, annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.25, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)
dev.off()




#------------- CD8 NP16 1131-TP-1 -------------------------

# separate TCR sequencing repeat sample 
# Peng carried out targeted scTCR sequencing -> bcl2fastq -> new MiXCR analyze amplicon
# filter by cloneCounts to get top alpha and beta chains per cell
# assess similarity between alpha and beta CDR3s 
# most similar CDR3s shared between alpha and beta in same cell = top alpha-beta pair for that cell

setwd('/t1-data/user/lfelce/TCR_analysis/1131-TP-1_cd8_new')

library(stringr)
library(data.table)
library(circlize)

# read in TRA names and convert to list
tra_names <- fread('/t1-data/user/lfelce/TCR_analysis/1131-TP-1_cd8_new/tra_names.txt', stringsAsFactors = F, header=F)

tra_names <- as.list(as.data.frame(t(tra_names)))

# remove ./ at start of name and replace with file path
tra_names <- tra_names %>% str_replace("./*", "")

# read in TRB names and convert to list
trb_names <- fread('/t1-data/user/lfelce/TCR_analysis/1131-TP-1_cd8_new/trb_names.txt', stringsAsFactors = F, header=F)

trb_names <- as.list(as.data.frame(t(trb_names)))

# remove ./ at start of name and replace with file path
trb_names <- trb_names %>% str_replace("./*", "")

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
  dat <- data.frame(c(mixcr_a[[i]][3],mixcr_a[[i]][4],mixcr_a[[i]][6],mixcr_a[[i]][7], mixcr_a[[i]][8]))
  dat$i <- i # keep track of which iteration produced it
  datalist[[i]] <- dat # add it to list
}
# combine columns for each cell, select chains with clone count more than 5
big_data = do.call(rbind, datalist)
tra <- big_data
tra <- tra[tra$Read.count > 5,]
colnames(tra) <- c("clone_count", "clone_fraction", "CDR3_aa", "TRAV", "TRAJ", "cell_number")
tra <- merge(tra, mixcr_a_names, by="cell_number")

# mixcr_b
datalist = list()
for (i in (1:length(mixcr_b))) {
  dat <- data.frame(c(mixcr_b[[i]][3],mixcr_b[[i]][4],mixcr_b[[i]][6],mixcr_b[[i]][7], mixcr_b[[i]][8]))
  dat$i <- i # keep track of which iteration produced it
  datalist[[i]] <- dat # add it to list
}
# combine columns for each cell, select chains with clone count more than 5
big_data = do.call(rbind, datalist)
trb <- big_data
trb <- trb[trb$Read.count > 5,]
colnames(trb) <- c("clone_count", "clone_fraction", "CDR3_aa", "TRBV", "TRBJ", "cell_number")
trb <- merge(trb, mixcr_b_names, by="cell_number")

# combine TRA and TRB dataframes
cd8_1131_TP1 <- merge(tra,trb, by="cell_name")

cd8_1131_TP1 <- mutate(cd8_1131_TP1, alpha=paste(TRAV, TRAJ, sep="_"))

cd8_1131_TP1 <- mutate(cd8_1131_TP1, beta=paste(TRBV, TRBJ, sep="_"))

library(stringdist)

# calculate similarity score using Levenshtein distance
shared_cdr3 <- stringsim(cd8_1131_TP1$CDR3_aa.x, cd8_1131_TP1$CDR3_aa.y, method="lv")

# add to data frame as separate column
cd8_1131_TP1$sim_score <- shared_cdr3                              

# create list of cell  names
list <- mixcr_a_names[,2]

# loop - select rows containing same cell name and make new object dat
# re-order rows by similarity score putting highest one at the top
# select top row only for each cell, make list of cells
datalist = list()

for (i in (1:length(list))) {
  dat <- cd8_1131_TP1[cd8_1131_TP1$cell_name == list[i],]
  dat <- dat[order(dat$sim_score, decreasing=TRUE),]
  dat <- dat[1,]
  datalist[[i]] <- dat # add it to list
}

# make one big data frame where each row is cell with highest similarity score for a-b pair
big_data = do.call(rbind, datalist)

# top alpha-beta pairing
top_ab <- as.data.frame(table(big_data$alpha, big_data$beta))
top_ab_mat <- as.matrix(as.data.frame.matrix(table(big_data$alpha, big_data$beta)))

setwd('/t1-data/user/lfelce/TCR_analysis/new_mixcr_results')
# circle plot
circos.clear()
set.seed(999)
pdf('cd8_1131-TP1_new_chorddiagram.pdf', width = 16, height = 12, useDingbats = FALSE)
chordDiagram(top_ab_mat, annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.25, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)
dev.off()


--------------------------------

# previous attempt - just selecting top alpha and top beta row per cell

setwd('/t1-data/user/lfelce/TCR_analysis/1131-TP-1_cd8_new')
listFiles = list.files()

# revise the structure of input files
for(i in 1:length(listFiles))
{
  DT = fread(listFiles[i])
  write.table(DT,listFiles[i],quote = F, row.names = F, sep = '\t')
}

library(stringr)

# read in TRA names and convert to list
tra_names <- fread('/t1-data/user/lfelce/TCR_analysis/1131-TP-1_cd8_new/tra_names.txt', stringsAsFactors = F, header=F)

tra_names <- as.list(as.data.frame(t(tra_names)))

# remove ./ at start of name
tra_names <- tra_names %>% str_replace("./*", "")

# read in TRB names and convert to list
trb_names <- fread('/t1-data/user/lfelce/TCR_analysis/1131-TP-1_cd8_new/trb_names.txt', stringsAsFactors = F, header=F)

trb_names <- as.list(as.data.frame(t(trb_names)))

# remove ./ at start of name
trb_names <- trb_names %>% str_replace("./*", "")


# parse mixcr files
# for some reason MiXCR makes clone count numeric rather than integer
# have modified code for parse.mixcr and saved in separated file
# call function parse.mixcr.ling
source('/t1-data/user/lfelce/TCR_analysis/parse_mixcr.R')

# TRA: generate list of lists (like mixcr parsed into R)
datalist = list()
for (i in (1:length(tra_names))) {
  dat <- parse.mixcr.ling(tra_names[i])
  datalist[[i]] <- dat # add it to list
}

# rename with sample name
names(datalist) <- tra_names
mixcr_a <- datalist

# TRB: generate list of lists (like mixcr parsed into R)
datalist = list()
for (i in (1:length(trb_names))) {
  dat <- parse.mixcr.ling(trb_names[i])
  datalist[[i]] <- dat # add it to list
}

# rename with sample name
names(datalist) <- trb_names
mixcr_b <- datalist

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

# remove file extensions from names
mixcr_a_names$cell_name <- mixcr_a_names$cell_name %>% str_replace("_L001.clonotypes.TRA.txt", "")
mixcr_b_names$cell_name <- mixcr_b_names$cell_name %>% str_replace("_L001.clonotypes.TRB.txt", "")

# convert mixcr lists to dataframe with just V.gene and J.gene info

# mixcr_a - select top 2 rows for each cell
datalist = list()
for (i in (1:length(mixcr_a))) {
  dat1 <- data.frame(c(mixcr_a[[i]][3],mixcr_a[[i]][7], mixcr_a[[i]][8]))
  dat <- dat1[1,]
  dat$i <- i # keep track of which iteration produced it
  datalist[[i]] <- dat # add it to list
}
# combine columns
big_data = do.call(rbind, datalist)
tra <- big_data
colnames(tra) <- c("Read.count", "TRAV", "TRAJ", "cell_number")
tra <- merge(tra, mixcr_a_names, by="cell_number")


# mixcr_b - select top row for each cell
datalist = list()
for (i in (1:length(mixcr_b))) {
  dat1 <- data.frame(c(mixcr_b[[i]][3],mixcr_b[[i]][7], mixcr_b[[i]][8]))
  dat <- dat1[1,]
  dat$i <- i # keep track of which iteration produced it
  datalist[[i]] <- dat # add it to list
}
# combine columns for each cell
big_data = do.call(rbind, datalist)
trb <- big_data 
colnames(trb) <- c("Read.count", "TRBV", "TRBJ", "cell_number")
trb <- merge(trb, mixcr_b_names, by="cell_number")

# combine TRA and TRB dataframes
cd8_1131 <- merge(tra,trb, by="cell_name")

cd8_1131 <- mutate(cd8_1131, alpha=paste(TRAV, TRAJ, sep="_"))

cd8_1131 <- mutate(cd8_1131, beta=paste(TRBV, TRBJ, sep="_"))


# tabulate to get dominant alpha-beta pairing

library(circlize)

cd8_np16_1131_TP1 <- as.matrix(as.data.frame.matrix(table(cd8_1131$alpha, cd8_1131$beta)))
cd8_np16_1131_TP1_table <- as.data.frame(table(cd8_1131$alpha, cd8_1131$beta))

setwd("/t1-data/user/lfelce/TCR_analysis/new_mixcr_results/")

list <- c("cd8_np16_1131_TP1")

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
