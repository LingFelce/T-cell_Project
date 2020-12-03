# single cell RNA-Seq analysis in R Studio

# will use R Studio on cbrg cluster

# module load R-cbrg/curret
# module load rstudio
# rstudio

# need following libraries: SingleCellExperiment, scater, ggplot2, pheatmap, ggfortify

library(pcaMethods)
library(SC3)
library(scater)
library(SingleCellExperiment)
library(pheatmap)
library(mclust)

# cleaning expression matrix 
# cell filtering - no UMI so will have to use scater to conduct PCA on set of QC metrics and use automatic outlier detection to identify potentially problematic cells
# spike ins? Have to email Neil Ashley again?
# gene filtering - genes must be filtered after cell filtering as some genes may only be detected in poor quality cells 

# plot PCA - before log transformation
tmp <- runPCA(
  umi[endog_genes, ],
  exprs_values = "counts"
)
plotPCA(
    tmp,
    colour_by = "batch",
    size_by = "total_features_by_counts",
    shape_by = "individual"
)

# plot PCA after log transformation
tmp <- runPCA(
  umi[endog_genes, ],
  exprs_values = "logcounts_raw"
)
plotPCA(
    tmp,
    colour_by = "batch",
    size_by = "total_features_by_counts",
    shape_by = "individual"
)

# plot PCA after QC
tmp <- runPCA(
  umi.qc[endog_genes, ],
  exprs_values = "logcounts_raw"
)
plotPCA(
    tmp,
    colour_by = "batch",
    size_by = "total_features_by_counts",
    shape_by = "individual"
)

# t-SNE
set.seed(123456)
tmp <- runTSNE(
    umi.qc[endog_genes, ],
    exprs_values = "logcounts_raw",
    perplexity = 130
)
plotTSNE(
    tmp,
    colour_by = "batch",
    size_by = "total_features_by_counts",
    shape_by = "individual"
)

# normalisation (reads)

# biological analysis - clustering


