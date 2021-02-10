###### Grouping single cells by CDR3 similarity ######
library(tcR)
library(tidyverse)
library(data.table)
library(stringr)
library(plyr)
library(dplyr)

#-----CD8 NP16------

# calculating similarity - CDR3 alpha
np16 <- read.csv("/t1-data/user/lfelce/TCR_analysis/shared_cdr3_cd8_np16.csv", header=TRUE)
np16_a <- np16[,c(2, 3, 8)]
colnames(np16_a) <- c("sequence", "TRAV", "name")
np16_a <- np16_a[!duplicated(np16_a),]

np16_clones <- read.csv("/t1-data/user/lfelce/TCR_analysis/cd8_np16_tcr_clones.csv", header=TRUE)
np16_clones_a <- np16_clones[,c(2, 7, 8)]
colnames(np16_clones_a) <- c("name", "sequence", "TRAV")

alpha <- rbind(np16_a, np16_clones_a)  

datalist = list()
for (i in (1:length(alpha$sequence))) {
  dat <- data.frame(agrep(pattern=alpha$sequence[i], x=alpha$sequence, 
                          max.distance=0.5,value=TRUE, fixed=TRUE))
  dat$i <- i # keep track of which iteration produced it
  datalist[[i]] <- dat # add it to list
}
big_data = do.call(rbind, datalist)
colnames(big_data) <- c("sequence", "group")
alpha_shared <- merge(big_data, alpha, by="sequence")
alpha_shared <- alpha_shared[!duplicated(alpha_shared$sequence, alpha_shared$name),]



alpha_shared <- merge(np16_a, big_data, by.x="CDR3a", by.y="sequence")
alpha_shared <- alpha_shared[!duplicated(alpha_shared$CDR3a),]

np16_clones <- read.csv("/t1-data/user/lfelce/TCR_analysis/cd8_np16_tcr_clones.csv", header=TRUE)
np16_clones_a <- np16_clones[,c(2, 7, 8)]

datalist = list()
for (i in (1:length(np16_clones_a$CloneName))) {
  dat <- data.frame(agrep(pattern=np16_clones_a$CDR3_aa.x[i], x=alpha_shared$CDR3a, 
                          max.distance=0.5,value=TRUE, fixed=TRUE))
  dat$i <- i # keep track of which iteration produced it
  dat$CloneName <- np16_clones_a$CloneName[i]
  datalist[[i]] <- dat # add it to list
}
big_data = do.call(rbind, datalist)
alpha_clones <- big_data
colnames(alpha_clones) <- c("sequence", "group", "clone_name")
alpha_clones <- alpha_clones[!duplicated(alpha_clones$sequence),]

# calculating similarity - CDR3 beta
np16_b <- np16[,c(1, 5, 6, 8, 9)]
np16_b <- np16_b[!duplicated(np16_b$CDR3b),]

datalist = list()
for (i in (1:length(np16_b$X))) {
  dat <- data.frame(agrep(pattern=np16_b$CDR3b[i], x=np16_b$CDR3b, 
                          max.distance=0.5,value=TRUE, fixed=TRUE))
  dat$i <- i # keep track of which iteration produced it
  datalist[[i]] <- dat # add it to list
}
big_data = do.call(rbind, datalist)
colnames(big_data) <- c("sequence", "group")
beta_shared <- merge(np16_b, big_data, by.x="CDR3b", by.y="sequence")
beta_shared <- beta_shared[!duplicated(beta_shared$CDR3b),]


#----CD8 ORF3a-28--------

# calculating similarity - CDR3 alpha
orf3a <- read.csv("/t1-data/user/lfelce/TCR_analysis/shared_cdr3_cd8_orf.csv", header=TRUE)
orf3a_a <- orf3a[,c(1:3, 8, 9)]
orf3a_a <- orf3a_a[!duplicated(orf3a_a$CDR3a),]

datalist = list()
for (i in (1:length(orf3a_a$X))) {
  dat <- data.frame(agrep(pattern=orf3a_a$CDR3a[i], x=orf3a_a$CDR3a, 
                          max.distance=0.5,value=TRUE, fixed=TRUE))
  dat$i <- i # keep track of which iteration produced it
  datalist[[i]] <- dat # add it to list
}
big_data = do.call(rbind, datalist)
colnames(big_data) <- c("sequence", "group")
alpha_shared <- merge(orf3a_a, big_data, by.x="CDR3a", by.y="sequence")
alpha_shared <- alpha_shared[!duplicated(alpha_shared$CDR3a),]

# calculating similarity - CDR3 beta
orf3a_b <- orf3a[,c(1, 5, 6, 8, 9)]
orf3a_b <- orf3a_b[!duplicated(orf3a_b$CDR3b),]

datalist = list()
for (i in (1:length(orf3a_b$X))) {
  dat <- data.frame(agrep(pattern=orf3a_b$CDR3b[i], x=orf3a_b$CDR3b, 
                          max.distance=0.5,value=TRUE, fixed=TRUE))
  dat$i <- i # keep track of which iteration produced it
  datalist[[i]] <- dat # add it to list
}
big_data = do.call(rbind, datalist)
colnames(big_data) <- c("sequence", "group")
beta_shared <- merge(orf3a_b, big_data, by.x="CDR3b", by.y="sequence")
beta_shared <- beta_shared[!duplicated(beta_shared$CDR3b),]
