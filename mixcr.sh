## below is script for mixcr. 

## create script folder in output directory

## copy and paste below into command line, it will generate individual script files for each sample

## then run second script at bottom, which will run all individual script files wherever they are saved.

cd /t1-data/user/lfelce/MiXCR/CD4_input2/ #input directory
DIR=/t1-data/user/lfelce/MiXCR/CD4_output/ # output directory

for NAME in $(find . -name '*_R1_001.fastq.gz' -printf "%f\n" | sed 's/_R1_001.fastq.gz//'); do # remove common ending of name
 
echo "$NAME"

p1='_R1_001.fastq.gz'

echo -e '#!/bin/sh
#$ -cwd
#$ -q batchq
module add mixcr
cd /t1-data/user/lfelce/MiXCR/CD4_input2/

mixcr align -p rna-seq -s hsa -OallowPartialAlignments=true' $NAME$p1 $DIR$NAME'.vdjca
# option preserves partial alignments for further use in assemblePartial

mixcr assemblePartial' $DIR$NAME'.vdjca' $DIR$NAME'_Rescued_it1.vdjca
mixcr assemblePartial' $DIR$NAME'_Rescued_it1.vdjca' $DIR$NAME'_Rescued_it2.vdjca
mixcr assemblePartial' $DIR$NAME'_Rescued_it2.vdjca' $DIR$NAME'_Rescued.vdjca
# several iterations to obtain more reads containing full CDR3 sequence

mixcr assemble -OaddReadsCountOnClustering=true -ObadQualityThreshold=15' $DIR$NAME'_Rescued.vdjca' $DIR$NAME'.clns
# decrease input quality threshold to 15 for poor quality data

mixcr exportClones' $DIR$NAME'.clns' $DIR$NAME'.txt' > $DIR'script/'$NAME'.sh'

done

#------------- run scripts on server - this should be separate script -------------------#

## copy and paste below into command line

cd /t1-data/user/lfelce/MiXCR/CD4_output/script/

for line in $(ls *.sh); do
sbatch $line
done

squeue -u lfelce


# create soft links for .txt files and move into folder for R Studio analysis

# find /t1-data/user/lfelce/MiXCR/CD4_output -name "*.txt" | xargs -I v_f ln -s v_f
# find /t1-data/user/lfelce/MiXCR/CD8_output -name "*.txt" | xargs -I v_f ln -s v_f

# check how many patient samples have valid clonotypes from mixcr output
# sort .txt files by size, delete all files that are 412 bytes (no clonotypes)
# copy and paste code below, it will print file names of remaining files in folder to console
# copy and paste this into a .txt file

for NAME in $(find . -name '*_R1_001.fastq.gz' -printf "%f\n" | sed 's/_R1_001.fastq.gz//'); do # remove common ending of name
 
echo "$NAME"

done

# count how many occurrences of each patient ID in sample file name to work out how many cells have clonotypes from each patient (out of 96)
grep -o -i 22 cd4_names.txt | wc -l
grep -o -i 25 cd4_names.txt | wc -l
grep -o -i 1062 cd4_names.txt | wc -l
grep -o -i 1493 cd4_names.txt | wc -l
grep -o -i 1504 cd4_names.txt | wc -l
grep -o -i 1525 cd4_names.txt | wc -l
