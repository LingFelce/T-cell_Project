##### Cell Ranger pipeline for analysing 10x data #####

## bcl conversion to fastq ##

cellranger mkfastq --help

Run Illumina demultiplexer on sample sheets that contain 10x-specific sample 
index sets, and generate 10x-specific quality metrics after the demultiplex.  
Any bcl2fastq argument works, except a few that are set by the pipeline 
to ensure proper trimming and sample indexing. The FASTQ output generated 
is the same as when running bcl2fastq directly.

# cellranger_mkfastq.sh

#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=mkfastq
#SBATCH --nodes=1
#SBATCH --mem=128G
#SBATCH --time=07-00:00:00
#SBATCH --output=%j_%x.out
#SBATCH --error=%j_%x.err


cd /t1-data/user/lfelce/10x_DONG171220/

module load cellranger/5.0.0
module load bcl2fastq

cellranger mkfastq --id=Dong_171220_10x \
--run=/t1-data/user/lfelce/10x_DONG171220/ \
--csv=/t1-data/user/lfelce/10x_DONG171220/Dong_171220_sample_sheet.csv

# cellranger_count.sh

cellranger count takes FASTQ files from cellranger mkfastq and performs alignment, filtering, barcode counting, and UMI counting. 
It uses the Chromium cellular barcodes to generate feature-barcode matrices, determine clusters, and perform gene expression analysis. 
The count pipeline can take input from multiple sequencing runs on the same GEM well. cellranger count also processes Feature Barcode data alongside Gene Expression reads.

#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=count
#SBATCH --nodes=1
#SBATCH --mem=128G
#SBATCH --time=07-00:00:00
#SBATCH --output=%j_%x.out
#SBATCH --error=%j_%x.err

cd /t1-data/user/lfelce/10x_DONG171220/Dong_171220_10x/outs/fastq_path/

module load cellranger/5.0.0

cellranger count --id=count_Dong_171220_10x --fastqs=/t1-data/user/lfelce/10x_DONG171220/Dong_171220_10x/outs/fastq_path/ --sample=Dong_171220_GEX --transcriptome=/databank/10x-rangers/refdata-cellranger-GRCh38-3.0.0

# cellranger aggr required, only for combining multiple runs of cellranger count to normalise runs to same sequencing depth.
# can use output from cellranger count filtered_feature_bc_matrix into Seurat

# For generating a hashtag count matrix from FASTQ files, please refer to https://github.com/Hoohm/CITE-seq-Count - luckily this has been installed on the cluster
# https://hoohm.github.io/CITE-seq-Count/Running-the-script/
-cbf 1 -cbl 16 -umif 17 -umil 26 # positions of the cellular and UMI barcodes for 10x

#### if running all 4 lanes together without merging, then do as 1 script ####
# cite-seq-count.sh

#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=cite-seq-count
#SBATCH --nodes=1
#SBATCH --mem=128G
#SBATCH --time=07-00:00:00
#SBATCH --output=%j_%x.out
#SBATCH --error=%j_%x.err

module load cite-seq-count

cd /t1-data/user/lfelce/10x_DONG171220/Dong_171220_10x/outs/fastq_path/

CITE-seq-Count -R1 Dong_171220_GEX_S1_L001_R1_001.fastq.gz,Dong_171220_GEX_S1_L002_R1_001.fastq.gz,Dong_171220_GEX_S1_L003_R1_001.fastq.gz,Dong_171220_GEX_S1_L004_R1_001.fastq.gz \
-R2 Dong_171220_GEX_S1_L001_R2_001.fastq.gz,Dong_171220_GEX_S1_L002_R2_001.fastq.gz,Dong_171220_GEX_S1_L003_R2_001.fastq.gz,Dong_171220_GEX_S1_L004_R2_001.fastq.gz \
-t hashtag_barcodes.csv \
-cbf 1 -cbl 16 -umif 17 -umil 26 -cells 2000 \
-o ./hashtag_count

###############

# cite-seq-count gave a very sparse matrix, very few hashtags mapped to cell barcodes/UMIs. 
# repeated using Feature Barcode extra options for cell ranger count
# libraries.csv - 1 row, specifying file path, sample name and library_type - Antibody Capture
# feature_ref - hashtag barcodes, TotalSeqC

# for some reason script doesn't work on CCB cluster - keeps saying no input fastq files - probably too many other fastq files in folder?
# copied GEX and Cell_Surface fastq files over to BMRC cluster and ran on there instead

# script below for combined Gene Expression and Antibody Capture

#!/bin/bash
#$ -wd /well/jknight/users/jln789/10x_Dong171220/
#$ -q short.qc

cd /well/jknight/users/jln789/10x_Dong171220/

module load CellRanger/5.0.0

cellranger count --id=combined_counts_2 \
--transcriptome=/well/htseq/Genomes/refdata-cellranger-3.1.0/GRCh38-combat/ \
--libraries=libraries.csv \
--feature-ref=feature_ref.csv


# script below for Antibody Capture only
#!/bin/bash
#$ -wd /well/jknight/users/jln789/10x_Dong171220/
#$ -q short.qc

cd /well/jknight/users/jln789/10x_Dong171220/

module load CellRanger/5.0.0

cellranger count --id=ab_counts_2 \
--transcriptome=/well/htseq/Genomes/refdata-cellranger-2020-A/refdata-gex-GRCh38-2020-A \
--libraries=libraries_2.csv \
--feature-ref=feature_ref_2.csv

# libraries.csv
fastqs,sample,library_type
/gpfs2/well/jknight/users/jln789/10x_Dong171220/,Dong_171220_Cell_Surface,Antibody Capture
/gpfs2/well/jknight/users/jln789/10x_Dong171220/,Dong_171220_GEX,Gene Expression


