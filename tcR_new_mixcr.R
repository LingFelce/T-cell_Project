##### Updated MiXCR-SmartSeq.sh output####
# Isar modified MiXCR script so that .txt output for each cell is TRA.txt, TRB.txt, TRD.txt, TRG.txt file
# read in TRA.txt and TRB.txt separately as mixcr list of dataframes
# pull out V.gene and J.gene info into 2 separate TRA and TRB dataframes
# combine information

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

#------------- MAIN ANALYSIS -------------

#------------- Read one MiXCR output
library(data.table)

#------------- Parse folder with MiXCR files.

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

# convert mixcr lists to dataframe with just V.gene and J.gene info

# mixcr_a
mixcr_a_names <- as.data.frame(names(mixcr_a))

datalist = list()
for (i in (1:180)) {
  dat <- data.frame(c(mixcr_a[[i]][7], mixcr_a[[i]][8]))
  dat$i <- i # maybe you want to keep track of which iteration produced it?
  datalist[[i]] <- dat # add it to your list
}
# combine columns for each cell, select only cells with only 1 row
big_data = do.call(rbind, datalist)
tra <- big_data %>% group_by(i) %>% filter(n() == 1)
colnames(tra) <- c("TRAV", "TRAJ", "A_cell_number")

# mixcr_b
mixcr_b_names <- as.data.frame(names(mixcr_b))

datalist = list()
for (i in (1:221)) {
  dat <- data.frame(c(mixcr_b[[i]][7], mixcr_b[[i]][8]))
  dat$i <- i # maybe you want to keep track of which iteration produced it?
  datalist[[i]] <- dat # add it to your list
}
# combine columns for each cell, select only cells with only 1 row
big_data = do.call(rbind, datalist)
trb <- big_data %>% group_by(i) %>% filter(n() == 1)
colnames(trb) <- c("TRBV", "TRBJ", "B_cell_number")


