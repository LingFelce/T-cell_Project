########## Isar's modified script for MiXCR ##############
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
# sorted files by size and removed empty output files 412 bytes - 590 files remaining
# need to read in .TRA.txt and .TRB.txt files separately

for NAME in $(find . -name '*TRA.txt'); do
echo "$NAME"
done

# copy and save output to cd8_np16_tra_names.txt
# repeat for TRB.txt
