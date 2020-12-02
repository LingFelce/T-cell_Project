## New MiXCR output shared CDR3 ##

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


library(data.table)

#------------------ CD8 NP16 ---------------------

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

# read in metadata
all_metadata = fread('/t1-data/user/lfelce/TCR_analysis/metadata_cd8.txt', stringsAsFactors = F)
colnames(all_metadata)

# change 5 to 005
for (i in (1:192)) {all_metadata[i,2] <- "005"}

# metadata has all samples - but not all samples have valid clones.

# select only file names which have valid clones
metadata_a <- all_metadata[is.element(all_metadata$sample, mixcr_a_names$cell_name),]
metadata_b <- all_metadata[is.element(all_metadata$sample, mixcr_b_names$cell_name),]

## alpha chain ##

# replace file names with patient index 
names(mixcr_a) = metadata_a$patient[match(names(mixcr_a),metadata_a$sample)]
length(unique(names(mixcr_a)))

mixcr_ordered = mixcr_a[order(names(mixcr_a))]
names(mixcr_ordered)

class(mixcr_a)

states = unique(names(mixcr_a))
o=1

mixcr_per_state = list()
for (o in 1:length(states)) {
  tem = mixcr_ordered[which(names(mixcr_ordered) %in% states[o])]
  tem2 <- do.call(rbind.data.frame, tem)
  dim(tem2)
  mixcr_per_state[[o]] = tem2
  length(mixcr_per_state)
  names(mixcr_per_state)[o] = states[o]
}
names(mixcr_per_state)
length(mixcr_per_state)

df <- do.call(rbind.data.frame, mixcr_per_state)
dim(df)

# shared CDR3 sequences

setwd('/t1-data/user/lfelce/TCR_analysis/new_mixcr_results')

imm.shared <- shared.repertoire(.data = mixcr_per_state, .type = 'avrc', .min.ppl = 2, .verbose = F)
head(imm.shared)

write.table(imm.shared, 'shared_repertoire_per_patient_alpha.txt', quote = F, row.names = F, sep = '\t')


## beta chain ##

# replace file names with patient index
names(mixcr_b) = metadata_b$patient[match(names(mixcr_b),metadata_b$sample)]
length(unique(names(mixcr_b)))

mixcr_ordered = mixcr_b[order(names(mixcr_b))]
names(mixcr_ordered)

class(mixcr_b)

states = unique(names(mixcr_b))
o=1

mixcr_per_state = list()
for (o in 1:length(states)) {
  tem = mixcr_ordered[which(names(mixcr_ordered) %in% states[o])]
  tem2 <- do.call(rbind.data.frame, tem)
  dim(tem2)
  mixcr_per_state[[o]] = tem2
  length(mixcr_per_state)
  names(mixcr_per_state)[o] = states[o]
}
names(mixcr_per_state)
length(mixcr_per_state)

df <- do.call(rbind.data.frame, mixcr_per_state)
dim(df)

# shared CDR3 sequences

setwd('/t1-data/user/lfelce/TCR_analysis/new_mixcr_results')

imm.shared <- shared.repertoire(.data = mixcr_per_state, .type = 'avrc', .min.ppl = 2, .verbose = F)
head(imm.shared)

write.table(imm.shared, 'shared_repertoire_per_patient_beta.txt', quote = F, row.names = F, sep = '\t')


#------------- CD 0RF3a-28 ------------------------

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

# read in metadata
all_metadata = fread('/t1-data/user/lfelce/TCR_analysis/metadata_cd8.txt', stringsAsFactors = F)
colnames(all_metadata)

# change 5 to 005
for (i in (1:192)) {all_metadata[i,2] <- "005"}

# metadata has all samples - but not all samples have valid clones.

# select only file names which have valid clones
metadata_a <- all_metadata[is.element(all_metadata$sample, mixcr_a_names$cell_name),]
metadata_b <- all_metadata[is.element(all_metadata$sample, mixcr_b_names$cell_name),]

## alpha chain ##

# replace file names with patient index
names(mixcr_a) = metadata_a$patient[match(names(mixcr_a),metadata_a$sample)]
length(unique(names(mixcr_a)))

mixcr_ordered = mixcr_a[order(names(mixcr_a))]
names(mixcr_ordered)

class(mixcr_a)

states = unique(names(mixcr_a))
o=1

mixcr_per_state = list()
for (o in 1:length(states)) {
  tem = mixcr_ordered[which(names(mixcr_ordered) %in% states[o])]
  tem2 <- do.call(rbind.data.frame, tem)
  dim(tem2)
  mixcr_per_state[[o]] = tem2
  length(mixcr_per_state)
  names(mixcr_per_state)[o] = states[o]
}
names(mixcr_per_state)
length(mixcr_per_state)

df <- do.call(rbind.data.frame, mixcr_per_state)
dim(df)

# shared CDR3 sequences

setwd('/t1-data/user/lfelce/TCR_analysis/new_mixcr_results')

imm.shared <- shared.repertoire(.data = mixcr_per_state, .type = 'avrc', .min.ppl = 2, .verbose = F)
head(imm.shared)

write.table(imm.shared, 'shared_repertoire_per_patient_alpha.txt', quote = F, row.names = F, sep = '\t')


## beta chain ##

# replace file names with patient index
names(mixcr_b) = metadata_b$patient[match(names(mixcr_b),metadata_b$sample)]
length(unique(names(mixcr_b)))

mixcr_ordered = mixcr_b[order(names(mixcr_b))]
names(mixcr_ordered)

class(mixcr_b)

states = unique(names(mixcr_b))
o=1

mixcr_per_state = list()
for (o in 1:length(states)) {
  tem = mixcr_ordered[which(names(mixcr_ordered) %in% states[o])]
  tem2 <- do.call(rbind.data.frame, tem)
  dim(tem2)
  mixcr_per_state[[o]] = tem2
  length(mixcr_per_state)
  names(mixcr_per_state)[o] = states[o]
}
names(mixcr_per_state)
length(mixcr_per_state)

df <- do.call(rbind.data.frame, mixcr_per_state)
dim(df)

# shared CDR3 sequences

setwd('/t1-data/user/lfelce/TCR_analysis/new_mixcr_results')

imm.shared <- shared.repertoire(.data = mixcr_per_state, .type = 'avrc', .min.ppl = 2, .verbose = F)
head(imm.shared)

write.table(imm.shared, 'shared_repertoire_per_patient_beta.txt', quote = F, row.names = F, sep = '\t')


#---------------- CD4 S34 + M24 ---------------

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

all_metadata = fread('/t1-data/user/lfelce/TCR_analysis/metadata_cd4.txt', stringsAsFactors = F)
colnames(all_metadata)
all_metadata$patient <- as.character(all_metadata$patient)

# rename patients 022 and 025
i=1
for (i in (1:96)) {all_metadata[i,2] <- "022"}
for (i in (193:288)) {all_metadata[i,2] <- "025"}

# metadata has all samples - but not all samples have valid clones.

# select only file names which have valid clones
metadata_a <- all_metadata[is.element(all_metadata$sample, mixcr_a_names$cell_name),]
metadata_b <- all_metadata[is.element(all_metadata$sample, mixcr_b_names$cell_name),]

## alpha chain ##

# replace file names with patient index
names(mixcr_a) = metadata_a$patient[match(names(mixcr_a),metadata_a$sample)]
length(unique(names(mixcr_a)))

mixcr_ordered = mixcr_a[order(names(mixcr_a))]
names(mixcr_ordered)

class(mixcr_a)

states = unique(names(mixcr_a))
o=1

mixcr_per_state = list()
for (o in 1:length(states)) {
  tem = mixcr_ordered[which(names(mixcr_ordered) %in% states[o])]
  tem2 <- do.call(rbind.data.frame, tem)
  dim(tem2)
  mixcr_per_state[[o]] = tem2
  length(mixcr_per_state)
  names(mixcr_per_state)[o] = states[o]
}
names(mixcr_per_state)
length(mixcr_per_state)

df <- do.call(rbind.data.frame, mixcr_per_state)
dim(df)

# shared CDR3 sequences

setwd('/t1-data/user/lfelce/TCR_analysis/new_mixcr_results')

imm.shared <- shared.repertoire(.data = mixcr_per_state, .type = 'avrc', .min.ppl = 2, .verbose = F)
head(imm.shared)

write.table(imm.shared, 'shared_repertoire_per_patient_alpha.txt', quote = F, row.names = F, sep = '\t')


## beta chain ##

# replace file names with patient index
names(mixcr_b) = metadata_b$patient[match(names(mixcr_b),metadata_b$sample)]
length(unique(names(mixcr_b)))

mixcr_ordered = mixcr_b[order(names(mixcr_b))]
names(mixcr_ordered)

class(mixcr_b)

states = unique(names(mixcr_b))
o=1

mixcr_per_state = list()
for (o in 1:length(states)) {
  tem = mixcr_ordered[which(names(mixcr_ordered) %in% states[o])]
  tem2 <- do.call(rbind.data.frame, tem)
  dim(tem2)
  mixcr_per_state[[o]] = tem2
  length(mixcr_per_state)
  names(mixcr_per_state)[o] = states[o]
}
names(mixcr_per_state)
length(mixcr_per_state)

df <- do.call(rbind.data.frame, mixcr_per_state)
dim(df)

# shared CDR3 sequences

setwd('/t1-data/user/lfelce/TCR_analysis/new_mixcr_results')

imm.shared <- shared.repertoire(.data = mixcr_per_state, .type = 'avrc', .min.ppl = 2, .verbose = F)
head(imm.shared)

write.table(imm.shared, 'shared_repertoire_per_patient_beta.txt', quote = F, row.names = F, sep = '\t')
