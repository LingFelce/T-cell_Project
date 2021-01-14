# Analysis of CD8 NP16 TCR clones from targeted TCR sequencing
# 50,000 cells per well, each well should be different clone
# Downloaded data from BaseSpace, converted to .fastq.gz
# Used MiXCR to align to TCR sequences - each well has separate TRA.txt and TRB.txt file
# Need to read in separate files - create list of dataframes for tra and trb

library(data.table)
library(tidyverse)

########## CD8 NP16 CLONES #################

# list of all .TRA.txt files
tra_list <- list.files(path = "./cd8_np16_tcr_clones", recursive = TRUE,
                       pattern = "\\TRA.txt$", 
                       full.names = TRUE)
trb_list <- list.files(path = "./cd8_np16_tcr_clones", recursive = TRUE,
                       pattern = "\\TRB.txt$", 
                       full.names = TRUE)

# Read all the files and create a FileName column to store filenames
tra_dt <- rbindlist(sapply(tra_list, fread, simplify = FALSE),
                    use.names = TRUE, idcol = "FileName")
trb_dt <- rbindlist(sapply(trb_list, fread, simplify = FALSE),
                    use.names = TRUE, idcol = "FileName")

# select columns FileName, cloneId, cloneCount, cloneFraction, aaSeqCDR3
tra <- tra_dt[,c("FileName", "cloneId", "cloneCount", "cloneFraction", "aaSeqCDR3")]
names(tra)[names(tra) == "aaSeqCDR3"] <- "CDR3_alpha_aa"

trb <- trb_dt[,c("FileName", "cloneId", "cloneCount", "cloneFraction", "aaSeqCDR3")]
names(trb)[names(trb) == "aaSeqCDR3"] <- "CDR3_beta_aa"

# filter based on clone count > 5
tra <- tra[tra$cloneCount > 5,]
trb <- trb[trb$cloneCount > 5,]

# import in clone names
cd8_np16_tra_clone_names <- read.csv("cd8_np16_tra_clone_names.csv")
cd8_np16_trb_clone_names <- read.csv("cd8_np16_trb_clone_names.csv")

# merge file names and clone names
tra <- merge(tra,cd8_np16_tra_clone_names, by="FileName")
trb <- merge(trb,cd8_np16_trb_clone_names, by="FileName")

# tidy up tables
tra <- tra[, c("CloneName", "cloneCount", "CDR3_alpha_aa")]
trb <- trb[, c("CloneName", "cloneCount", "CDR3_beta_aa")]

# probably easier to export files separately as tra has 411 rows and trb has 216 rows (most dual alpha single beta?)
write.csv(tra, "cd8_np16_tra_cdr3.csv")
write.csv(trb, "cd8_np16_trb_cdr3.csv")

# some CDR3 sequences are shared across different patients - tabulate and export
alpha_freq <- as.data.frame(table(tra$CDR3_alpha_aa))
shared_cdr3_alpha <- merge(tra, alpha_freq, by.x="CDR3_alpha_aa", by.y="Var1")

beta_freq <- as.data.frame(table(trb$CDR3_beta_aa))
shared_cdr3_beta <- merge(trb, beta_freq, by.x="CDR3_beta_aa", by.y="Var1")

write.csv(shared_cdr3_alpha, "cd8_np16_shared_alpha_cdr3_tcr.csv")
write.csv(shared_cdr3_beta, "cd8_np16_shared_beta_cdr3_tcr.csv")

#------- Repeat using tcR to get CDR3, V and J info ------

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

setwd('/t1-data/user/lfelce/TCR_analysis/cd8_np16_tcr_clones')

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

setwd('/t1-data/user/lfelce/TCR_analysis')
# create list of file name and file number - will be same for both (1-96)
clone_names <- read.csv("cd8_np16_tra_clone_names.csv")
clone_names <-tibble::rownames_to_column(clone_names, "cell_number")

# mixcr a
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
tra <- merge(clone_names, tra, by="cell_number")

# mixcr b
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
trb <- merge(clone_names, trb, by="cell_number")

# calculate clone_fraction as percentage
tra$clone_percentage <- tra$clone_fraction * 100
tra <- tra %>% filter(clone_percentage > 1)

trb$clone_percentage <- trb$clone_fraction * 100
trb <- trb %>% filter(clone_percentage > 1)

write.csv(tra, "cd8_np16_tra_cdr3.csv")
write.csv(trb, "cd8_np16_trb_cdr3.csv")

# merge V and J columns and create new alpha/beta column
tra <- mutate(tra, alpha=paste(TRAV, TRAJ, sep="_"))
trb <- mutate(trb, beta=paste(TRBV, TRBJ, sep="_"))

# merge alpha and beta
cd8_np16 <- merge(tra,trb, by="CloneName")

write.csv(cd8_np16, "cd8_np16_tcr_clones.csv")

# some CDR3 sequences are shared across different patients - tabulate and export
alpha_freq <- as.data.frame(table(tra$CDR3_aa))
shared_cdr3_alpha <- merge(tra, alpha_freq, by.x="CDR3_aa", by.y="Var1")

beta_freq <- as.data.frame(table(trb$CDR3_aa))
shared_cdr3_beta <- merge(trb, beta_freq, by.x="CDR3_aa", by.y="Var1")

write.csv(shared_cdr3_alpha, "cd8_np16_shared_alpha_cdr3_tcr.csv")
write.csv(shared_cdr3_beta, "cd8_np16_shared_beta_cdr3_tcr.csv")

# from merged table
tcr_freq <- as.data.frame(table(cd8_np16$alpha, cd8_np16$beta))
shared_tcr <- merge(cd8_np16, tcr_freq, by.x=c("alpha","beta"), by.y=c("Var1","Var2"))

write.csv(shared_tcr, "cd8_np16_clones_shared_tcr.csv")

################### CD8 ORF3a-28 TCR CLONES ############################

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

# Nucleotide sequence and CDR3 position of each gene segment
genesegments$TRBV[1:10,]

setwd('/t1-data/user/lfelce/TCR_analysis/cd8_orf3a-28_tcr_clones')

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

setwd('/t1-data/user/lfelce/TCR_analysis')
# create list of file name and file number - will be same for both (1-96)
clone_names <- read.csv("cd8_orf3a-28_clone_names.csv")
clone_names <-tibble::rownames_to_column(clone_names, "cell_number")

# mixcr a
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
tra <- merge(clone_names, tra, by="cell_number")

# mixcr b
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
trb <- merge(clone_names, trb, by="cell_number")

tra$clone_percentage <- tra$clone_fraction * 100
tra <- tra %>% filter(clone_percentage > 1)

trb$clone_percentage <- trb$clone_fraction * 100
trb <- trb %>% filter(clone_percentage > 1)

write.csv(tra, "cd8_orf3a-28_tra_cdr3.csv")
write.csv(trb, "cd8_orf3a-28_trb_cdr3.csv")

# merge V and J columns and create new alpha/beta column
tra <- mutate(tra, alpha=paste(TRAV, TRAJ, sep="_"))
trb <- mutate(trb, beta=paste(TRBV, TRBJ, sep="_"))

# merge alpha and beta
cd8_orf3a <- merge(tra,trb, by="CloneName")

write.csv(cd8_orf3a, "cd8_orf3a-28_tcr_clones.csv")

# some CDR3 sequences are shared across different patients - tabulate and export
# alpha_freq <- as.data.frame(table(tra$CDR3_aa))
# shared_cdr3_alpha <- merge(tra, alpha_freq, by.x="CDR3_aa", by.y="Var1")
# 
# beta_freq <- as.data.frame(table(trb$CDR3_aa))
# shared_cdr3_beta <- merge(trb, beta_freq, by.x="CDR3_aa", by.y="Var1")
# 
# write.csv(shared_cdr3_alpha, "cd8_orf3a-28_shared_alpha_cdr3_tcr.csv")
# write.csv(shared_cdr3_beta, "cd8_orf3a-28_shared_beta_cdr3_tcr.csv")

# from merged table
tcr_freq <- as.data.frame(table(cd8_orf3a$alpha, cd8_orf3a$beta))
shared_tcr <- merge(cd8_orf3a, tcr_freq, by.x=c("alpha","beta"), by.y=c("Var1","Var2"))

write.csv(shared_tcr, "cd8_orf3a-28_clones_shared_tcr.csv")
