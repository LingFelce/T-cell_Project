##### Updated MiXCR-SmartSeq.sh output####
# Isar modified MiXCR script so that .txt output for each cell is TRA.txt, TRB.txt, TRD.txt, TRG.txt file
library(tcR)
library(tidyverse)
library(data.table)
library(stringr)
library(plyr)
library(dplyr)


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


#------------- CD8 NP16 -------------------------

setwd('/t1-data/user/lfelce/TCR_analysis/cd8_np16_new/')

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
# removed bulk samples for now
# keep all cells

# mixcr_a
datalist = list()
for (i in (1:length(mixcr_a))) {
  dat <- data.frame(c(mixcr_a[[i]][3], mixcr_a[[i]][4], mixcr_a[[i]][6],mixcr_a[[i]][7], mixcr_a[[i]][8]))
  dat$i <- i # keep track of which iteration produced it
  datalist[[i]] <- dat # add it to list
}
# combine columns for each cell and filter to keep cells with 1 or 2 alphas
big_data = do.call(rbind, datalist)
tra <- big_data 
colnames(tra) <- c("clone_count", "clone_fraction", "CDR3_alpha", "TRAV", "TRAJ", "cell_number")
tra <- merge(tra, mixcr_a_names, by="cell_number")


# mixcr_b
datalist = list()
for (i in (1:length(mixcr_b))) {
  dat <- data.frame(c(mixcr_b[[i]][3], mixcr_b[[i]][4], mixcr_b[[i]][6], mixcr_b[[i]][7], mixcr_b[[i]][8]))
  dat$i <- i # keep track of which iteration produced it
  datalist[[i]] <- dat # add it to list
}
# combine columns for each cell
big_data = do.call(rbind, datalist)
trb <- big_data 
colnames(trb) <- c("clone_count", "clone_fraction","CDR3_beta", "TRBV", "TRBJ", "cell_number")
trb <- merge(trb, mixcr_b_names, by="cell_number")

# export alpha and beta lists separately
# write.csv(tra, "/t1-data/user/lfelce/TCR_analysis/new_mixcr_results/all_alpha_cd8_np16.csv")
# write.csv(trb, "/t1-data/user/lfelce/TCR_analysis/new_mixcr_results/all__beta_cd8_np16.csv")

np16_a <- tra
np16_b <- trb


#------------- CD8 ORF3a-28 -------------------------

setwd('/t1-data/user/lfelce/TCR_analysis/cd8_orf_new/')

tra_list <- list.files(path = ".", recursive = TRUE,
                       pattern = "\\TRA.txt$", 
                       full.names = TRUE)
tra_list <- tra_list %>% str_replace("./*", "")
tra_list <-tra_list[!str_detect(tra_list,pattern="minibulk")]

trb_list <- list.files(path = ".", recursive = TRUE,
                       pattern = "\\TRB.txt$", 
                       full.names = TRUE)
trb_list <- trb_list %>% str_replace("./*", "")
trb_list <-trb_list[!str_detect(trb_list,pattern="minibulk")]

# parse mixcr files
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

# mixcr_a
datalist = list()
for (i in (1:length(mixcr_a))) {
  dat <- data.frame(c(mixcr_a[[i]][3], mixcr_a[[i]][4], mixcr_a[[i]][6],mixcr_a[[i]][7], mixcr_a[[i]][8]))
  dat$i <- i # keep track of which iteration produced it
  datalist[[i]] <- dat # add it to list
}
# combine columns for each cell and filter to keep cells with 1 or 2 alphas
big_data = do.call(rbind, datalist)
tra <- big_data 
colnames(tra) <- c("clone_count", "clone_fraction", "CDR3_alpha", "TRAV", "TRAJ", "cell_number")
tra <- merge(tra, mixcr_a_names, by="cell_number")


# mixcr_b
datalist = list()
for (i in (1:length(mixcr_b))) {
  dat <- data.frame(c(mixcr_b[[i]][3], mixcr_b[[i]][4], mixcr_b[[i]][6], mixcr_b[[i]][7], mixcr_b[[i]][8]))
  dat$i <- i # keep track of which iteration produced it
  datalist[[i]] <- dat # add it to list
}
# combine columns for each cell
big_data = do.call(rbind, datalist)
trb <- big_data 
colnames(trb) <- c("clone_count", "clone_fraction","CDR3_beta", "TRBV", "TRBJ", "cell_number")
trb <- merge(trb, mixcr_b_names, by="cell_number")

# export alpha and beta lists separately
# write.csv(tra, "/t1-data/user/lfelce/TCR_analysis/new_mixcr_results/all_alpha_cd8_orf3a-28.csv")
# write.csv(trb, "/t1-data/user/lfelce/TCR_analysis/new_mixcr_results/all_beta_cd8_orf3a-28.csv")

orf3a_a <- tra
orf3a_b <- trb

#------------- CD4 S34 and M24 -------------------------

setwd('/t1-data/user/lfelce/TCR_analysis/cd4_new/')
tra_list <- list.files(path = ".", recursive = TRUE,
                       pattern = "\\TRA.txt$", 
                       full.names = TRUE)
tra_list <- tra_list %>% str_replace("./*", "")
tra_list <-tra_list[!str_detect(tra_list,pattern="minibulk")]

trb_list <- list.files(path = ".", recursive = TRUE,
                       pattern = "\\TRB.txt$", 
                       full.names = TRUE)
trb_list <- trb_list %>% str_replace("./*", "")
trb_list <-trb_list[!str_detect(trb_list,pattern="minibulk")]

# parse mixcr files
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

# mixcr_a
datalist = list()
for (i in (1:length(mixcr_a))) {
  dat <- data.frame(c(mixcr_a[[i]][3], mixcr_a[[i]][4], mixcr_a[[i]][6],mixcr_a[[i]][7], mixcr_a[[i]][8]))
  dat$i <- i # keep track of which iteration produced it
  datalist[[i]] <- dat # add it to list
}
# combine columns for each cell and filter to keep cells with 1 or 2 alphas
big_data = do.call(rbind, datalist)
tra <- big_data 
colnames(tra) <- c("clone_count", "clone_fraction", "CDR3_alpha", "TRAV", "TRAJ", "cell_number")
tra <- merge(tra, mixcr_a_names, by="cell_number")


# mixcr_b
datalist = list()
for (i in (1:length(mixcr_b))) {
  dat <- data.frame(c(mixcr_b[[i]][3], mixcr_b[[i]][4], mixcr_b[[i]][6], mixcr_b[[i]][7], mixcr_b[[i]][8]))
  dat$i <- i # keep track of which iteration produced it
  datalist[[i]] <- dat # add it to list
}
# combine columns for each cell
big_data = do.call(rbind, datalist)
trb <- big_data 
colnames(trb) <- c("clone_count", "clone_fraction","CDR3_beta", "TRBV", "TRBJ", "cell_number")
trb <- merge(trb, mixcr_b_names, by="cell_number")

# # export alpha and beta lists separately
# write.csv(tra, "/t1-data/user/lfelce/TCR_analysis/new_mixcr_results/all_alpha_cd4_s34_m24.csv")
# write.csv(trb, "/t1-data/user/lfelce/TCR_analysis/new_mixcr_results/all_beta_cd4_s34_m24.csv")

cd4_a <- tra
cd4_b <- trb

#----------------MERGE------------------

alpha <- rbind(np16_a, orf3a_a, cd4_a)
beta <- rbind(np16_b, orf3a_b, cd4_b)

alpha <- alpha[,-c(2,3)]
beta <- beta[,-c(2,3)]

alpha$patient <- ifelse(grepl("005", alpha$cell_name), "1",
                 ifelse(grepl("1131-TP-1", alpha$cell_name), "Dong231120",
                 ifelse(grepl("1131-TP-2", alpha$cell_name),"2",
                 ifelse(grepl("1153", alpha$cell_name), "3", 
                 ifelse(grepl("1201", alpha$cell_name), "4",
                 ifelse(grepl("1105", alpha$cell_name), "1",
                 ifelse(grepl("1134", alpha$cell_name), "2", 
                 ifelse(grepl("1525-TP-1", alpha$cell_name), "Dong121020", 
                 ifelse(grepl("1525-TP-2", alpha$cell_name), "3", 
                 ifelse(grepl("022", alpha$cell_name), "1", 
                 ifelse(grepl("1062", alpha$cell_name),"2", 
                 ifelse(grepl("025", alpha$cell_name), "1", 
                 ifelse(grepl("1493", alpha$cell_name), "2", 
                 ifelse(grepl("1504", alpha$cell_name), "3", 
                 ifelse(grepl("1525_M24", alpha$cell_name), "4", "")))))))))))))))

alpha$epitope <- ifelse(grepl("NP16|B7_SPR", alpha$cell_name), "NP16",
                        ifelse(grepl("ORF3a-28", alpha$cell_name), "ORF3a-28",
                               ifelse(grepl("M24", alpha$cell_name), "M24",
                                      ifelse(grepl("S34", alpha$cell_name), "S34", ""))))

alpha$subtype <- ifelse(grepl("NP16|ORF3a-28", alpha$epitope), "CD8",
                        ifelse(grepl("M24|S34", alpha$epitope), "CD4", ""))

beta$patient <- ifelse(grepl("005", beta$cell_name), "1",
                ifelse(grepl("1131-TP-1", beta$cell_name), "Dong231120",
                ifelse(grepl("1131-TP-2", beta$cell_name),"2",
                ifelse(grepl("1153", beta$cell_name), "3", 
                ifelse(grepl("1201", beta$cell_name), "4",
                ifelse(grepl("1105", beta$cell_name), "1",
                ifelse(grepl("1134", beta$cell_name), "2", 
                ifelse(grepl("1525-TP-1", beta$cell_name), "Dong121020", 
                ifelse(grepl("1525-TP-2", beta$cell_name), "3", 
                ifelse(grepl("022", beta$cell_name), "1", 
                ifelse(grepl("1062", beta$cell_name),"2", 
                ifelse(grepl("025", beta$cell_name), "1", 
                ifelse(grepl("1493", beta$cell_name), "2", 
                ifelse(grepl("1504", beta$cell_name), "3", 
                ifelse(grepl("1525_M24", beta$cell_name), "4", "")))))))))))))))

beta$epitope <- ifelse(grepl("NP16|B7_SPR", beta$cell_name), "NP16",
                        ifelse(grepl("ORF3a-28", beta$cell_name), "ORF3a-28",
                               ifelse(grepl("M24", beta$cell_name), "M24",
                                      ifelse(grepl("S34", beta$cell_name), "S34", ""))))

beta$subtype <- ifelse(grepl("NP16|ORF3a-28", beta$epitope), "CD8",
                        ifelse(grepl("M24|S34", beta$epitope), "CD4", ""))

write.csv(alpha, "/t1-data/user/lfelce/TCR_analysis/new_mixcr_results/all_alpha.csv")
write.csv(beta, "/t1-data/user/lfelce/TCR_analysis/new_mixcr_results/all_beta.csv")
