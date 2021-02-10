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


datalist = list()
for (i in (1:length(np16))) {
  dat <- data.frame(agrep(pattern=np16$CDR3a[i], x=np16$CDR3a, 
                          max.distance=0.3,value=TRUE, fixed=TRUE))
  dat$i <- i # keep track of which iteration produced it
  datalist[[i]] <- dat # add it to list
}
big_data = do.call(rbind, datalist)
colnames(big_data) <- c("sequence", "group")
alpha_shared <- merge(np16, big_data, by.x="CDR3a", by.y="sequence")


datalist = list()
for (i in (1:length(np16))) {
  dat <- data.frame(agrep(pattern=np16$CDR3b[i], x=np16$CDR3b, 
                          max.distance=0.3,value=TRUE, fixed=TRUE))
  dat$i <- i # keep track of which iteration produced it
  datalist[[i]] <- dat # add it to list
}
big_data = do.call(rbind, datalist)
colnames(big_data) <- c("sequence", "group")
beta_shared <- merge(np16, big_data, by.x="CDR3b", by.y="sequence")
