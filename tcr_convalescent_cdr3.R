# Comparing 6 month convalescent CDR3 sequences to acute/1 month convalescent CDR3 sequences

library(tidyverse)
library(stringr)
library(data.table)

###### NP16 beta #####

# smartseq <- read.csv("/t1-data/user/lfelce/TCR_analysis/new_mixcr_results/cd8_np16_complete_sc_tcr.csv")
# smartseq_b <- smartseq[,c("cell_name", "CDR3_beta", "beta")]
# smartseq_b <- smartseq_b[!duplicated(smartseq_b$cell_name),]
# smartseq_b$CDR3_beta <- paste("TRB:", smartseq_b$CDR3_beta, sep="")
# 
# sc10x <- read.csv("/t1-data/user/lfelce/CellRanger_VDJ/10x_Dong050121_TCR/TCR_T2/outs/clonotypes.csv")
# sc10x_b <- sc10x[sc10x$cdr3s_aa %like% "TRB", ]
# 
# sc10x_b2 <- as.data.frame(str_split_fixed(sc10x_b$cdr3s_aa, ";", 4))
# sc10x_b3 <- cbind(sc10x_b, sc10x_b2)
# 
# combine <- merge(sc10x_b3, smartseq_b, by.x="V1", by.y="CDR3_beta")
# combine <- combine[,c(1:5, 12, 13)]
# colnames(combine) <- c("shared_CDR3_beta", "clonotype_id", "frequency", "proportion", 
#                        "clonotype_cdr3", "cell_name", "cell_trb")
# combine <- combine[,c(6, 7, 1:5)]
# 
# combine2 <- merge(sc10x_b3, smartseq_b, by.x="V2", by.y="CDR3_beta")
# combine2 <- combine2[,c(1:5, 12, 13)]
# colnames(combine2) <- c("shared_CDR3_beta", "clonotype_id", "frequency", "proportion", 
#                         "clonotype_cdr3", "cell_name", "cell_trb")
# combine2 <- combine2[,c(6, 7, 1:5)]
# 
# np16_beta <- rbind(combine, combine2)
# 
# setwd("/t1-data/user/lfelce/TCR_analysis/shared_convalescent_cdr3s/")
# 
# write.csv(np16_beta, "np16_beta_shared_convalescent_cdr3s.csv")
# 
# ###### NP16 alpha #####
# 
# smartseq <- read.csv("/t1-data/user/lfelce/TCR_analysis/new_mixcr_results/cd8_np16_complete_sc_tcr.csv")
# smartseq_a <- smartseq[,c("cell_name", "CDR3_alpha", "alpha")]
# smartseq_a <- smartseq_a[!duplicated(smartseq_a$cell_name),]
# smartseq_a$CDR3_alpha <- paste("TRA:", smartseq_a$CDR3_alpha, sep="")
# 
# sc10x <- read.csv("/t1-data/user/lfelce/CellRanger_VDJ/10x_Dong050121_TCR/TCR_T2/outs/clonotypes.csv")
# sc10x_a <- sc10x[sc10x$cdr3s_aa %like% "TRA", ]
# 
# sc10x_a2 <- as.data.frame(str_split_fixed(sc10x_b$cdr3s_aa, ";", 4))
# sc10x_a3 <- cbind(sc10x_a, sc10x_a2)
# 
# combine <- merge(sc10x_a3, smartseq_a, by.x="V2", by.y="CDR3_alpha")
# combine <- combine[,c(1:5, 12, 13)]
# colnames(combine) <- c("shared_CDR3_alpha", "clonotype_id", "frequency", "proportion", 
#                        "clonotype_cdr3", "cell_name", "cell_tra")
# combine <- combine[,c(6, 7, 1:5)]
# 
# combine2 <- merge(sc10x_a3, smartseq_a, by.x="V3", by.y="CDR3_alpha")
# combine2 <- combine2[,c(1:5, 12, 13)]
# colnames(combine2) <- c("shared_CDR3_alpha", "clonotype_id", "frequency", "proportion", 
#                         "clonotype_cdr3", "cell_name", "cell_tra")
# combine2 <- combine2[,c(6, 7, 1:5)]
# 
# combine3 <- merge(sc10x_a3, smartseq_a, by.x="V4", by.y="CDR3_alpha")
# combine3 <- combine3[,c(1:5, 12, 13)]
# colnames(combine3) <- c("shared_CDR3_alpha", "clonotype_id", "frequency", "proportion", 
#                         "clonotype_cdr3", "cell_name", "cell_tra")
# combine3 <- combine3[,c(6, 7, 1:5)]
# 
# np16_alpha <- rbind(combine, combine2, combine3)
# 
# write.csv(np16_alpha, "np16_alpha_shared_convalescent_cdr3s.csv")

## ORF3a-28 beta ##

smartseq <- read.csv("/t1-data/user/lfelce/TCR_analysis/new_mixcr_results/cd8_orf3a-28_complete_sc_tcr.csv")
smartseq_b <- smartseq[,c("cell_name", "CDR3_beta", "beta")]
smartseq_b <- smartseq_b[!duplicated(smartseq_b$cell_name),]
smartseq_b$CDR3_beta <- paste("TRB:", smartseq_b$CDR3_beta, sep="")

sc10x <- read.csv("/t1-data/user/lfelce/CellRanger_VDJ/10x_Dong050121_TCR/TCR_T2/outs/clonotypes.csv")
sc10x_b <- sc10x[sc10x$cdr3s_aa %like% "TRB", ]

sc10x_b2 <- as.data.frame(str_split_fixed(sc10x_b$cdr3s_aa, ";", 4))
sc10x_b3 <- cbind(sc10x_b, sc10x_b2)

# match <- as.data.frame(agrep(smartseq_b$CDR3_beta[1], sc10x_b3$V1, max=0.2, value=TRUE))

# compare each ORF3a-28 CDR3 beta sequence with 10x sequence
# max.distance 0.25 so each cell gets a match
# value = FALSE - gives row names from 10x dataframe, correlates with clonotype number

datalist = list()
for (i in (1:length(smartseq_b$cell_name))) {
  dat <- data.frame(agrep(pattern=smartseq_b$CDR3_beta[i], x=sc10x_b3$V1, 
                          max.distance=0.25,value=FALSE, fixed=TRUE))
  dat$i <- i # keep track of which iteration produced it
  datalist[[i]] <- dat # add it to list
}
big_data = do.call(rbind, datalist)
colnames(big_data) <- c("clonotype", "sc_number")
V1 <- big_data

datalist = list()
for (i in (1:length(smartseq_b$cell_name))) {
  dat <- data.frame(agrep(pattern=smartseq_b$CDR3_beta[i], x=sc10x_b3$V2, 
                          max.distance=0.25,value=FALSE, fixed=TRUE))
  dat$i <- i # keep track of which iteration produced it
  datalist[[i]] <- dat # add it to list
}
big_data = do.call(rbind, datalist)
colnames(big_data) <- c("clonotype", "sc_number")
V2 <- big_data

match <- rbind(V1, V2)

# Look for only unique clonotypes
# not matching each cell to each clonotype, will only select at clonotype level
match <- as.data.frame(unique(match$clonotype))



