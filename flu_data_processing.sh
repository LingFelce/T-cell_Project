###### FEATURE COUNTS FOR DONOR 2 #########

#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=featurecounts
#SBATCH --nodes=1
#SBATCH --mem=128G
#SBATCH --time=07-00:00:00
#SBATCH --output=%j_%x.out
#SBATCH --error=%j_%x.err

module load subread/2.0.0 

cd /t1-data/user/ypeng/P170665/fastq_files/

featureCounts -p -a /databank/igenomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf -o /t1-data/user/lfelce/210121_counts.txt ./*.bam

########## FEATURE COUNTS FOR DONOR 1 ##############

module load subread/2.0.0 

cd /t1-data/user/ypeng/P170335/downloaded_files/lane2/

featureCounts -p -a /databank/igenomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf -o /t1-data/user/lfelce/D1_210121_counts.txt ./*.bam




########## OLD CODE, NOT SURE IF CORRECT FILES! ###############

# soft link all .fastq.gz files from Peng's folder (lane 1 and lane 2) to my folder
find /t1-data/user/ypeng/P170335/downloaded_files/lane1 -name "*.fastq.gz" | xargs -I v_f ln -s v_f
find /t1-data/user/ypeng/P170335/downloaded_files/lane2 -name "*.fastq.gz" | xargs -I v_f ln -s v_f

# find unique identifier/name for each sample and save in file called ID (should have 192 IDs per lane)
ls -1 *_1*.gz | awk -F '_' '{print $3}' | sort | uniq > ID

# merge files from lane 1 (WTCHG_390445) and lane 2 (WTCHG_392409) into a new folder called merge and save file name after unique ID. 
# do separately for read 1 (_1) and read 2 (_2)
for i in `cat ./ID`; do cat WTCHG_390445_$i\_1.fastq.gz WTCHG_392409_$i\_1.fastq.gz > merge/$i\_1.fastq.gz; done
for i in `cat ./ID`; do cat WTCHG_390445_$i\_2.fastq.gz WTCHG_392409_$i\_2.fastq.gz > merge/$i\_2.fastq.gz; done
