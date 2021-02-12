###### Grouping single cells by CDR3 similarity ######
library(tcR)
library(tidyverse)
library(data.table)
library(stringr)
library(plyr)
library(dplyr)

#-----CD8 NP16------

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
# keep cells with 1 or 2 alphas
# keep cells with 1 beta
# merge alpha and beta and keep cells that express only 1 alpha or only 1 beta

# alpha

# mixcr_a
datalist = list()
for (i in (1:length(mixcr_a))) {
  dat <- data.frame(c(mixcr_a[[i]][3], mixcr_a[[i]][4], mixcr_a[[i]][6],mixcr_a[[i]][7], mixcr_a[[i]][8]))
  dat$i <- i # keep track of which iteration produced it
  datalist[[i]] <- dat # add it to list
}
# combine columns for each cell and filter to keep cells with 1 or 2 alphas
big_data = do.call(rbind, datalist)
tra <- big_data %>% group_by(i) %>% filter(n() <= 2)
colnames(tra) <- c("clone_count", "clone_fraction", "CDR3_alpha", "TRAV", "TRAJ", "cell_number")
tra <- merge(tra, mixcr_a_names, by="cell_number")
tra <- mutate(tra, alpha=paste(TRAV, TRAJ, sep="_"))
tra <- tra[,c("cell_name", "CDR3_alpha", "alpha")]
colnames(tra) <- c("name", "CDR3a", "TRA")
tra <- unique(tra)

# clones
np16_clones <- read.csv("/t1-data/user/lfelce/TCR_analysis/cd8_np16_tcr_clones.csv", header=TRUE)
np16_clones_a <- np16_clones[,c("CloneName", "CDR3_aa.x", "alpha")]
colnames(np16_clones_a) <- c("name", "CDR3a", "TRA")
np16_clones_a <- unique(np16_clones_a)

alpha <- rbind(tra, np16_clones_a)  

datalist = list()
for (i in (1:length(alpha$name))) {
  dat <- data.frame(agrep(pattern=alpha$CDR3a[i], x=alpha$CDR3a, 
                          max.distance=0.3,value=TRUE, fixed=TRUE))
  dat$i <- i # keep track of which iteration produced it
  datalist[[i]] <- dat # add it to list
}
big_data = do.call(rbind, datalist)
colnames(big_data) <- c("CDR3a", "group")
alpha_shared <- merge(big_data, alpha, by="CDR3a")
alpha_shared <- alpha_shared[order(alpha_shared$group),]
alpha_shared <- alpha_shared[!duplicated(alpha_shared$name),]
alpha_shared <- alpha_shared %>% group_by(group) %>% filter(n() >= 2)

write.csv(alpha_shared, "/t1-data/user/lfelce/TCR_analysis/cdr3_similarity/cd8_np16_sc_clone_shared_cdr3.csv")

# beta

# mixcr_b
datalist = list()
for (i in (1:length(mixcr_b))) {
  dat <- data.frame(c(mixcr_b[[i]][3], mixcr_b[[i]][4], mixcr_b[[i]][6], mixcr_b[[i]][7], mixcr_b[[i]][8]))
  dat$i <- i # keep track of which iteration produced it
  datalist[[i]] <- dat # add it to list
}
# combine columns for each cell
big_data = do.call(rbind, datalist)
trb <- big_data %>% group_by(i) %>% filter(n() == 1)
colnames(trb) <- c("clone_count", "clone_fraction","CDR3_beta", "TRBV", "TRBJ", "cell_number")
trb <- merge(trb, mixcr_b_names, by="cell_number")
trb <- mutate(trb, beta=paste(TRBV, TRBJ, sep="_"))
trb <- trb[,c("cell_name", "CDR3_beta", "beta")]
colnames(trb) <- c("name", "CDR3b", "TRB")
trb <- unique(trb)

np16_clones_b <- np16_clones[,c("CloneName", "CDR3_aa.y", "beta")]
colnames(np16_clones_b) <- c("name", "CDR3b", "TRB")
np16_clones_b <- unique(np16_clones_b)

beta <- rbind(trb, np16_clones_b)  

datalist = list()
for (i in (1:length(beta$name))) {
  dat <- data.frame(agrep(pattern=beta$CDR3b[i], x=beta$CDR3b, 
                          max.distance=0.3,value=TRUE, fixed=TRUE))
  dat$i <- i # keep track of which iteration produced it
  datalist[[i]] <- dat # add it to list
}
big_data = do.call(rbind, datalist)
colnames(big_data) <- c("CDR3b", "group")
beta_shared <- merge(big_data, beta, by="CDR3b")
beta_shared <- beta_shared[order(beta_shared$group),]
beta_shared <- beta_shared[!duplicated(beta_shared$name),]
beta_shared <- beta_shared %>% group_by(group) %>% filter(n() >= 2)

write.csv(beta_shared, "/t1-data/user/lfelce/TCR_analysis/cdr3_similarity/cd8_np16_sc_clone_shared_cdr3.csv")





#----CD8 ORF3a-28--------

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

# beta only

# mixcr_b
datalist = list()
for (i in (1:length(mixcr_b))) {
  dat <- data.frame(c(mixcr_b[[i]][3], mixcr_b[[i]][4], mixcr_b[[i]][6], mixcr_b[[i]][7], mixcr_b[[i]][8]))
  dat$i <- i # keep track of which iteration produced it
  datalist[[i]] <- dat # add it to list
}
# combine columns for each cell
big_data = do.call(rbind, datalist)
trb <- big_data %>% group_by(i) %>% filter(n() == 1)
colnames(trb) <- c("clone_count", "clone_fraction","CDR3_beta", "TRBV", "TRBJ", "cell_number")
trb <- merge(trb, mixcr_b_names, by="cell_number")
trb <- mutate(trb, beta=paste(TRBV, TRBJ, sep="_"))
trb <- trb[,c("cell_name", "CDR3_beta", "beta")]
colnames(trb) <- c("name", "CDR3b", "TRB")
trb <- unique(trb)

# clones
orf3a_clones <- read.csv("/t1-data/user/lfelce/TCR_analysis/cd8_orf3a-28_tcr_clones.csv", header=TRUE)
orf3a_clones <- orf3a_clones %>% filter(clone_percentage.y > 10)
orf3a_clones_b <- orf3a_clones[,c("CloneName", "CDR3_aa.y", "beta")]
colnames(orf3a_clones_b) <- c("name", "CDR3b", "TRB")
orf3a_clones_b <- unique(orf3a_clones_b)

beta <- rbind(trb, orf3a_clones_b)  

datalist = list()
for (i in (1:length(beta$name))) {
  dat <- data.frame(agrep(pattern=beta$CDR3b[i], x=beta$CDR3b, 
                          max.distance=0.3,value=TRUE, fixed=TRUE))
  dat$i <- i # keep track of which iteration produced it
  datalist[[i]] <- dat # add it to list
}
big_data = do.call(rbind, datalist)
colnames(big_data) <- c("CDR3b", "group")
beta_shared <- merge(big_data, beta, by="CDR3b")
beta_shared <- beta_shared[order(beta_shared$group),]
beta_shared <- beta_shared[!duplicated(beta_shared$name),]
beta_shared <- beta_shared %>% group_by(group) %>% filter(n() >= 2)

write.csv(beta_shared, "/t1-data/user/lfelce/TCR_analysis/cdr3_similarity/cd8_orf3a_sc_clone_shared_cdr3.csv")
