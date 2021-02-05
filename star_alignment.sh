# for CCB cluster #

#!/bin/sh
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --mem=128G
#SBATCH --time=07-00:00:00
module load STAR/2.6.1d
cd /t1-data/user/lfelce/STAR/MPE_TCR_clone/

STAR --runThreadN 10 --readFilesCommand gunzip \
-c --outSAMtype BAM SortedByCoordinate \
--genomeDir /databank/indices/star/hg19 \
--readFilesIn ./H09_S93_L001_R1_001.fastq.gz ./H09_S93_L001_R2_001.fastq.gz \
--outFileNamePrefix ./MRN345_Cl07_

# for BMRC cluster #

#!/bin/bash
#$ -wd /well/jknight/users/jln789/MPE_TCR_clone/
#$ -q short.qc
###$ -o output.log
###$ -e error.log

module load STAR/2.7.1a-foss-2018b                                                                                                                                                                                                           cd /well/jknight/users/jln789/MPE_TCR_clone/

STAR --runThreadN 10 --readFilesCommand gunzip \
-c --outSAMtype BAM SortedByCoordinate \
--genomeDir /well/htseq/Genomes/refdata-cellranger-2020-A/refdata-gex-GRCh38-2020-A/star \
--readFilesIn ./H09_S93_L001_R1_001.fastq.gz ./H09_S93_L001_R2_001.fastq.gz \
--outFileNamePrefix ./MRN345_Cl07_
