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

cellranger count --id=count_Dong_171220_10x \
--fastqs=/t1-data/user/lfelce/10x_DONG171220/Dong_171220_10x/outs/fastq_path/ \
--sample=Dong171220 \
--transcriptome=/databank/10x-rangers/refdata-cellranger-GRCh38-3.0.0

