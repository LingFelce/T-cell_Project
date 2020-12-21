########## Isar's modified script for MiXCR ##############

# see bottom of page for updated version with 2 reads and new batch queue information
# note: specific example below won't work as files are from TCR sequencing, not RNA-sequencing!

# copy and paste below into command line. 
# then when have generated scripts, copy and paste second bit into command line

cd /t1-data/user/lfelce/MiXCR/CD8_input/
DIR=/t1-data/user/lfelce/MiXCR/CD8_output_new/

for NAME in $(find . -name '*_R1_001.fastq.gz' -printf "%f\n" | sed 's/_R1_001.fastq.gz//'); do
 
echo "$NAME"

p1='_R1_001.fastq.gz' 

echo -e '#!/bin/sh
#$ -cwd
#$ -q batchq
module load mixcr
cd /t1-data/user/lfelce/MiXCR/CD8_output_new/
mixcr analyze shotgun -s hsa --starting-material rna --contig-assembly --only-productive /t1-data/user/lfelce/MiXCR/CD8_input/'$NAME$p1 $NAME '--receptor-type tcr' > $DIR'script/'$NAME'.sh'

done

#------------- run scripts on server
cd /t1-data/user/lfelce/MiXCR/CD8_output_new/script/

for line in $(ls *.sh); do
sbatch $line
done

squeue -u lfelce


################### preparing files for R Studio analysis ###################

# soft link to /t1-data/user/lfelce/TCR_analysis/cd8_np16_new/
find /t1-data/user/lfelce/MiXCR/CD8_output_new -name "*.txt" | xargs -I v_f ln -s v_f

# sorted files by size and removed empty output files 412 bytes - 590 files remaining
# need to read in .TRA.txt and .TRB.txt files separately

for NAME in $(find . -name '*TRA.txt'); do
echo "$NAME"
done

for NAME in $(find . -name '*TRB.txt'); do
echo "$NAME"
done

# copy and save output to cd8_np16_tra_names.txt
# repeat for TRB.txt

############# updated version of new MiXCR-SmartSeq.sh script ############
cd /t1-data/user/lfelce/MiXCR/1131-TP-1_CD8_new_input/
DIR=/t1-data/user/lfelce/MiXCR/1131-TP-1_CD8_new_output/

for NAME in $(find . -name '*_R1_001.fastq.gz' -printf "%f\n" | sed 's/_R1_001.fastq.gz//'); do

echo "$NAME"

p1='_R1_001.fastq.gz'
p2='_R2_001.fastq.gz'

echo -e '#! /bin/sh
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --mem=128G

module load mixcr

cd /t1-data/user/lfelce/MiXCR/1131-TP-1_CD8_new_output/

mixcr analyze shotgun -s hsa --starting-material rna --contig-assembly 
--only-productive /t1-data/user/lfelce/MiXCR/1131-TP-1_CD8_new_input/'$NAME$p1 
'/t1-data/user/lfelce/MiXCR/1131-TP-1_CD8_new_input/'$NAME$p2
$NAME '--receptor-type tcr' > 
