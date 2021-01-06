# Analysis of CD8 NP16 TCR clones from targeted TCR sequencing
# 50,000 cells per well, each well should be different clone
# Downloaded data from BaseSpace, converted to .fastq.gz
# Used MiXCR 3.0.3 to align to TCR sequences - each well has separate TRA.txt and TRB.txt file
# Need to read in separate files - create list of dataframes for tra and trb


library(data.table)

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
cd8_np16_tra_clone_names <- read.csv("/gpfs3/well/jknight/users/jln789/tcr_analysis/cd8_np16_tra_clone_names.csv")
cd8_np16_trb_clone_names <- read.csv("/gpfs3/well/jknight/users/jln789/tcr_analysis/cd8_np16_trb_clone_names.csv")

# merge file names and clone names
tra <- merge(tra,cd8_np16_tra_clone_names, by="FileName")
trb <- merge(trb,cd8_np16_trb_clone_names, by="FileName")
