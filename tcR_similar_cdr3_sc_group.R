###### Grouping single cells by CDR3 similarity ######
library(tcR)
library(tidyverse)
library(data.table)
library(stringr)
library(plyr)
library(dplyr)

## CD8 NP16 ##

# calculating similarity - CDR3 alpha
np16 <- read.csv("/t1-data/user/lfelce/TCR_analysis/shared_cdr3_cd8_np16.csv", header=TRUE)
np16_a <- np16[,c(1:3, 8, 9)]
np16_a <- np16_a[!duplicated(np16_a$CDR3a),]

datalist = list()
for (i in (1:length(np16_a$X))) {
  dat <- data.frame(agrep(pattern=np16_a$CDR3a[i], x=np16_a$CDR3a, 
                          max.distance=0.1,value=TRUE, fixed=TRUE))
  dat$i <- i # keep track of which iteration produced it
  datalist[[i]] <- dat # add it to list
}
big_data = do.call(rbind, datalist)
colnames(big_data) <- c("sequence", "group")
alpha_shared <- merge(np16_a, big_data, by.x="CDR3a", by.y="sequence")


np16_b <- np16[,c(1, 5, 6, 8, 9)]
np16_b <- np16_b[!duplicated(np16_b$CDR3b),]

datalist = list()
for (i in (1:length(np16_b$X))) {
  dat <- data.frame(agrep(pattern=np16_b$CDR3b[i], x=np16_b$CDR3b, 
                          max.distance=0.1,value=TRUE, fixed=TRUE))
  dat$i <- i # keep track of which iteration produced it
  datalist[[i]] <- dat # add it to list
}
big_data = do.call(rbind, datalist)
colnames(big_data) <- c("sequence", "group")
beta_shared <- merge(np16_b, big_data, by.x="CDR3b", by.y="sequence")
