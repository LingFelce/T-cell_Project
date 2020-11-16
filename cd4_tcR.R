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

setwd('/t1-data/user/lfelce/TCR_analysis/cd4/')
listFiles = list.files()

# revise the structure of input files
for(i in 1:length(listFiles))
{
  DT = fread(listFiles[i])
  write.table(DT,listFiles[i],quote = F, row.names = F, sep = '\t')
}

#read files
mixcr <- parse.folder("/t1-data/user/lfelce/TCR_analysis/cd4/", 'mixcr')

all_metadata = fread('/t1-data/user/lfelce/TCR_analysis/metadata_cd4.txt', stringsAsFactors = F)
colnames(all_metadata)

# metadata has all samples - but not all samples have valid clones.

# names of files in folder- 398 files
name_list <- fread('/t1-data/user/lfelce/TCR_analysis/cd4_names.txt', stringsAsFactors = F, header=F)

# select only file names which have valid clones
metadata <- all_metadata[is.element(all_metadata$sample, name_list$V1),]


# mixcr_backup = mixcr

# replace file names with patient index - should have 4 different patients
names(mixcr) = metadata$patient[match(names(mixcr),metadata$sample)]
length(unique(names(mixcr)))

mixcr_ordered = mixcr[order(names(mixcr))]
names(mixcr_ordered)

class(mixcr)

states = unique(names(mixcr))
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


#------------- proportion of the most abundant clonotypes'
top.proportion(mixcr_per_state, 10)

vis.top.proportions(mixcr_per_state)

setwd('/t1-data/user/lfelce/TCR_analysis/cd4_results/')
pdf('top_proportions_per_state.pdf', width = 12, height = 8, useDingbats = FALSE)
vis.top.proportions(mixcr_per_state)
dev.off()

#------------- Clonal space homeostasis - identification of hyper-expanded clonotypes
#Clonal space homeostasis is a useful statistic of how much space is occupied by clonotypes with specific proportions.

mixcr.space <- clonal.space.homeostasis(mixcr_per_state)
write.table(mixcr.space, 'mixcr_clonal_space_homeostasis_per_state.txt', quote = F, row.names = T, sep = '\t')

vis.clonal.space(mixcr.space)


pdf('clonal.space_per_state.pdf', width = 12, height = 8, useDingbats = FALSE)
vis.clonal.space(mixcr.space)
dev.off()


pdf('clonal.space_high_low_per_state.pdf', width = 12, height = 8, useDingbats = FALSE)
vis.clonal.space(clonal.space.homeostasis(mixcr_per_state, c(Low = .05, High = 1)))
dev.off()


#------------- Search for a target CDR3 sequences - CMV
# cmv <- data.frame(CDR3.amino.acid.sequence = c('ASSFDTEAF'), stringsAsFactors = F)
# 
# cmv
# 
# mixcr <- set.rank(mixcr)
# cmv.imm.ex <-
#   find.clonotypes(.data = mixcr, .targets = cmv[,1], .method = 'exact',
#                   .col.name = c('Read.count', 'Total.insertions'),
#                   .verbose = F)
# head(cmv.imm.ex)

#------------- Plot gene usage - select whole section as part of for bit
condition = names(mixcr_per_state)
i=1
for(i in condition)
{
  
  mixcr_subset = mixcr_per_state[grep(i, names(mixcr_per_state))]
  
  #---- Top cross
  # mixcr.top <- top.cross(.data = mixcr_subset, .n = seq(500, 10000, 500), .verbose = T, .norm = T)
  # 
  # pdf(paste0(i,'_mixcr_Topcross_per_state.pdf'), width = 25, height = 10, useDingbats = FALSE)
  # print(top.cross.plot(mixcr.top))
  # dev.off()
   
  mixcr.jusage <- geneUsage(mixcr_subset, HUMAN_TRAV)
  
  pdf(paste0(i,'_HUMAN_TRAV_mixcrJ-usage_dodge_per_state.pdf'), width = 12, height = 8, useDingbats = FALSE)
  print(vis.gene.usage(mixcr.jusage, .main = 'HUMAN_TRAV-usage', .dodge = T))
  dev.off()
  
  mixcr.jusage <- geneUsage(mixcr_subset, HUMAN_TRAJ)
  pdf(paste0(i,'_HUMAN_TRAJ_mixcrJ-usage_dodge_per_state.pdf'), width = 12, height = 8, useDingbats = FALSE)
  print(vis.gene.usage(mixcr.jusage, .main = 'HUMAN_TRAJ-usage', .dodge = T))
  dev.off()
  
  mixcr.jusage <- geneUsage(mixcr_subset, HUMAN_TRBJ)
  pdf(paste0(i,'_HUMAN_TRBJ_mixcrJ-usage_dodge_per_state.pdf'), width = 12, height = 8, useDingbats = FALSE)
  print(vis.gene.usage(mixcr.jusage, .main = 'HUMAN_TRBJ-usage', .dodge = T))
  dev.off()
  
  mixcr.vusage <- geneUsage(mixcr_subset, HUMAN_TRBV)
  pdf(paste0(i,'_HUMAN_TRBV_mixcrV-usage_dodge_per_state.pdf'), width = 30, height = 10, useDingbats = FALSE)
  print(vis.gene.usage(mixcr.vusage, .main = 'HUMAN_TRBV-usage', .dodge = T))
  dev.off()
  
  ### Logo-like plot
  #Logo-like graphs for visualisation of nucleotide or amino acid motif sequences / profiles.
  # km <- get.kmers(mixcr[[j]]$CDR3.amino.acid.sequence, .head = 100, .k = 7, .verbose = F)
  # d <- kmer.profile(km)
  # 
  # pdf(paste0(names(mixcr_subset)[j],'_mixcr_logo_motif_sequences_per_state.pdf'), width = 12, height = 8, useDingbats = FALSE)
  # print(vis.logo(d))
  # dev.off()
  
}

#--- Gene usage comparing
#Shannon entropy and Jensen-Shannon divergence
entropy(0:100, .laplace = 1)
imm.js <- js.div.seg(mixcr_per_state, HUMAN_TRBV, .verbose = F) 
pdf(paste0(i,'_mixcr_HUMAN_TRBV-usage_compare_states.pdf'), width = 30, height = 20, useDingbats = FALSE)
print(vis.radarlike(imm.js))
dev.off()

imm.js <- js.div.seg(mixcr_per_state, HUMAN_TRBJ, .verbose = F) 
pdf(paste0(i,'_mixcr_HUMAN_TRBJ-usage_compare_states.pdf'), width = 30, height = 20, useDingbats = FALSE)
print(vis.radarlike(imm.js))
dev.off()

imm.js <- js.div.seg(mixcr_per_state, HUMAN_TRAV, .verbose = F) 
pdf(paste0(i,'_mixcr_HUMAN_TRAV-usage_compare_states.pdf'), width = 30, height = 20, useDingbats = FALSE)
print(vis.radarlike(imm.js))
dev.off()

imm.js <- js.div.seg(mixcr_per_state, HUMAN_TRAJ, .verbose = F) 
pdf(paste0(i,'_mixcr_HUMAN_TRAJ-usage_compare_states.pdf'), width = 30, height = 20, useDingbats = FALSE)
print(vis.radarlike(imm.js))
dev.off()

imm.js <- js.div.seg(mixcr_per_state, HUMAN_IGHJ, .verbose = F) 
pdf(paste0(i,'_mixcr_HUMAN_IGHJ-usage_compare_states.pdf'), width = 30, height = 20, useDingbats = FALSE)
print(vis.radarlike(imm.js))
dev.off()

imm.js <- js.div.seg(mixcr_per_state, HUMAN_IGHV, .verbose = F) 
pdf(paste0(i,'_mixcr_HUMAN_IGHV-usage_compare_states.pdf'), width = 30, height = 20, useDingbats = FALSE)
print(vis.radarlike(imm.js))
dev.off()

#------------- Gene usage comparing PER CONDITION
#Compute entropy of V-segment usage for each data frame.
# entropy.seg(mixcr, HUMAN_TRBV)
# js.div.seg(mixcr[1:2], HUMAN_TRBV, .verbose = F)

imm.shared <- shared.repertoire(.data = mixcr_per_state, .type = 'avrc', .min.ppl = 2, .verbose = F)
head(imm.shared)
count=shared.representation(imm.shared)

write.table(imm.shared, 'SharedRepertoire_percondition_per_state.txt', quote = F, row.names = F, sep = '\t')
write.table(count, 'shared_representation_percondition_per_state.txt', quote = F, row.names = F, sep = '\t')


#------------- Principal Component Analysis (PCA)

pdf(paste0('mixcr_CD4_PCA_HUMAN_TRAV.pdf'), width = 12, height = 8, useDingbats = FALSE)
pca.segments(mixcr_per_state, .genes = HUMAN_TRAV, .text = T)  # Plot PCA results of V-segment usage.
dev.off()

pdf(paste0('mixcr_CD4_PCA_HUMAN_TRAJ.pdf'), width = 12, height = 8, useDingbats = FALSE)
pca.segments(mixcr_per_state, .genes = HUMAN_TRAJ, .text = T)  # Plot PCA results of V-segment usage.
dev.off()

pdf(paste0('mixcr_CD4_PCA_HUMAN_TRBV.pdf'), width = 12, height = 8, useDingbats = FALSE)
pca.segments(mixcr_per_state, .genes = HUMAN_TRBV, .text = T)  # Plot PCA results of V-segment usage.
dev.off()

pdf(paste0('mixcr_CD4_PCA_HUMAN_TRBJ.pdf'), width = 12, height = 8, useDingbats = FALSE)
pca.segments(mixcr_per_state, .genes = HUMAN_TRBJ, .text = T)  # Plot PCA results of V-segment usage.
dev.off()

pdf(paste0('mixcr_CD4_PCA_HUMAN_IGHJ.pdf'), width = 12, height = 8, useDingbats = FALSE)
pca.segments(mixcr_per_state, .genes = HUMAN_IGHJ, .text = T)  # Plot PCA results of V-segment usage.
dev.off()

pdf(paste0('mixcr_PCA_CD4_HUMAN_IGHV.pdf'), width = 12, height = 8, useDingbats = FALSE)
pca.segments(mixcr_per_state, .genes = HUMAN_IGHV, .text = T)  # Plot PCA results of V-segment usage.
dev.off()

#------------- Repertoire overlap analysis

# Equivalent to intersectClonesets(mixcr, "n0e", .norm = T)
# repOverlap(mixcr, 'exact', 'nuc', .norm = T, .verbose = F)
# Intersect by amino acid clonotypes + V genes
# repOverlap(mixcr, 'exact', 'aa', .vgene = T, .verbose = F)
# Plot a heatmap of the number of shared clonotypes.


pdf('heatmap_number_shared_clonotypes_exact_CDR3.nucleotide.sequence.pdf', width = 15, height = 10, useDingbats = FALSE)
print(vis.heatmap(repOverlap(mixcr_per_state, 'exact', 'nuc', .vgene = F, .verbose = T, .norm = T), .title = 'mixcr - (ave)-intersection', .labs = 'Samples', .signif.digits = 2, .size.text = 5, .no.labs = F, .scientific = T))
dev.off()

overlap = repOverlap(mixcr_per_state, 'exact', 'nuc', .vgene = T, .verbose = F, .norm = T)
write.table(overlap, 'repOverlap_exact_CDR3.nucleotide.sequence.txt', quote = F, row.names = F, sep = '\t')

pdf('heatmap_number_shared_clonotypes_exact_CDR3.amino_acids.sequence.pdf', width = 15, height = 10, useDingbats = FALSE)
print(vis.heatmap(repOverlap(mixcr_per_state, 'exact', 'aa', .vgene = F, .verbose = T, .norm = T), .title = 'mixcr - (ave)-intersection', .labs = 'Samples', .signif.digits = 2, .size.text = 5, .no.labs = F, .scientific = T))
dev.off()

overlap = repOverlap(mixcr_per_state, 'exact', 'aa', .vgene = T, .verbose = F, .norm = T)
write.table(overlap, 'repOverlap_exact_CDR3.amino_acids.sequence.txt', quote = F, row.names = F, sep = '\t')

pdf('heatmap_number_shared_clonotypes_exact_CDR3_jaccard.pdf', width = 15, height = 10, useDingbats = FALSE)
print(vis.heatmap(repOverlap(mixcr_per_state, 'jaccard', 'aa', 'read.count', .vgene = F, .verbose = T), .title = 'mixcr - (ave)-intersection', .labs = 'Samples', .signif.digits = 2, .size.text = 5, .no.labs = F, .scientific = T))
dev.off()
overlap = repOverlap(mixcr_per_state, 'jaccard', 'aa', 'read.count', .vgene = T, .verbose = F, .norm = T)
write.table(overlap, 'repOverlap_jaccard_CDR3.amino_acids.sequence.txt', quote = F, row.names = T, sep = '\t')


# overlap_sub = overlap[-56,]

# library(superheat)
# superheat(overlap, 
#           
#           # place dendrograms on columns and rows 
#           row.dendrogram = T, 
#           col.dendrogram = T,
#           
#           # make gridlines white for enhanced prettiness
#           grid.hline.col = "white",
#           grid.vline.col = "white",
#           
#           # rotate bottom label text
#           bottom.label.text.angle = 90,
#           legend.text.size = 15, left.label.text.size = 5, bottom.label.text.size =5) 
# 
# 
# View(overlap[1:10,1:10])
#------------- Overlap statistics and tests
#### Overlap Z-score (OZ-score) - a measure for ???abnormality??? in overlaps
# `ozScore`

#### Monte Carlo permutation test for pairwise and one-vs-all-wise within- and inter-group differences in a set of repertoires
# `permutDistTest`
# `pca2euclid`

### Shared repertoire
#To investigate a shared clonotype among several repertoires ("shared repertoire") the package provided the `shared.repertoire` function along with functions for computing the shared repertoire statistics. 
#The `shared.representation` function computes the number of shared clonotypes for each repertoire for each degree of sharing (i.e., number of people, in which indicated amount of clones have been found). 
#The function `shared.summary` is equivalent to `repOverlap(, 'exact')` but applies to the shared repertoire data frame. 
#Measuring distances among repertoires using the cosine similarity on vector of counts of shared sequences is also possible with the `cosine.sharing` function.
# Compute shared repertoire of amino acid CDR3 sequences and V genes
# which has been found in two or more people and return the Read.count column
# of such clonotypes from each data frame in the input list.
imm.shared <- shared.repertoire(.data = mixcr, .type = 'avrc', .min.ppl = 2, .verbose = F)
head(imm.shared)
shared.representation(imm.shared)  # Number of shared sequences.

write.table(imm.shared, 'SharedRepertoire_per_state.txt', quote = F, row.names = F, sep = '\t')


## Diversity evaluation
#For assessing the distribution of clonotypes in the given repertoire, *tcR* provides functions for evaluating the diversity (functions `diversity` and `inverse.simpson`) and the skewness of the clonal distribution (functions `gini` and `gini.simpson`), and a general interface to all of this functions `repDiversity`, which user should use to estimate the diversity of clonesets. 
#Function `diversity` (`repDiversity(your_clonesets, "div")`) computes the ecological diversity index (with parameter `.q` for penalties for clones with large count). 
#Function `inverse.simpson` (`repDiversity(your_clonesets, "inv.simp")`) computes the Inverse Simpson Index (i.e., inverse probability of choosing two similar clonotypes). 
#Function `gini` (`repDiversity(your_clonesets, "gini")`) computes the economical Gini index of clonal distribution. Function `gini.simpson` (`repDiversity(your_clonesets, "gini.simp")`) computes the Gini-Simpson index. 
#Function `chao1` (`repDiversity(your_clonesets, "chao1")`) computes the Chao1 index, its SD and two 95 perc CI. Function `repDiversity` accepts single clonesets as well as a list of clonesets. Parameter `.quant` specifies which column to use for computing the diversity (print `?repDiversity` to see more information about input arguments).
# Evaluate the diversity of clones by the ecological diversity index.
repDiversity(mixcr, 'div', 'read.count')
sapply(mixcr, function (x) diversity(x$Read.count))

# Compute the diversity as the inverse probability of choosing two similar clonotypes.
repDiversity(mixcr, 'inv.simp', 'read.prop')
sapply(mixcr, function (x) inverse.simpson(x$Read.proportion))

# Evaluate the skewness of clonal distribution.
repDiversity(mixcr, 'gini.simp', 'read.prop')
sapply(mixcr, function (x) gini.simpson(x$Read.proportion))

# Compute diversity of repertoire using Chao index.
repDiversity(mixcr, 'chao1', 'read.count')
sapply(mixcr, function (x) chao1(x$Read.count))

# output of diversity evaluation?

#------------- Mutation networks
#Mutation network (or a mutation graph) is a graph with vertices representing nucleotide or in-frame amino acid sequences (out-of-frame amino acid sequences will be automatically filtered out by *tcR* functions for mutation network creating) and edges which connecting pairs of sequences with hamming distance (parameter *.method* = 'hamm') or edit distance (parameter *.method* = 'lev') between them no more than specified in the *.max.errors* function parameter of the `mutation.network` function. To create a mutation network first what you need is to make a shared repertoires and then apply the `mutation.network` function to this shared repertoire:

## End
