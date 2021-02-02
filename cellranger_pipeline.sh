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

cd /t1-data/user/lfelce/10x_DONG171220/Dong_171220_10x/outs/fastq_path/
DIR=/t1-data/user/lfelce/10x_DONG171220/Dong_171220_10x/outs/fastq_path/hashtag_count/

for NAME in $(find . -name '*GEX_*_R1_001.fastq.gz' -printf "%f\n" | sed 's/_R1_001.fastq.gz//'); do

echo "$NAME"

p1='_R1_001.fastq.gz'
p2='_R2_001.fastq.gz'

echo -e '#! /bin/sh
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --mem=128G

module load cite-seq-count

cd /t1-data/user/lfelce/10x_DONG171220/Dong_171220_10x/outs/fastq_path/hashtag_count/

CITE-seq-Count -R1 /t1-data/user/lfelce/10x_DONG171220/Dong_171220_10x/outs/fastq_path/'$NAME$p1' -R2 /t1-data/user/lfelce/10x_DONG171220/Dong_171220_10x/outs/fastq_path/'$NAME$p2' -t hashtag_barcodes.csv -cbf 1 -cbl 16 -umif 17 -umil 26 -cells 2000 -o .' > $DIR'script/'$NAME'.sh'

done

#-----------------------

cd /t1-data/user/lfelce/10x_DONG171220/Dong_171220_10x/outs/fastq_path/hashtag_count/script/

for line in $(ls *.sh); do
sbatch $line
done

squeue -u lfelce

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
