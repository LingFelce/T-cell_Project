## below is script for mixcr. 

## create script folder in output directory

## copy and paste below into command line, it will generate individual script files for each sample

## then run second script at bottom, which will run all individual script files wherever they are saved.

cd /t1-data/user/lfelce/MiXCR/CD8_input/ #input directory
DIR=/t1-data/user/lfelce/MiXCR/CD8_output/ # output directory

for NAME in $(find . -name '*_R1_001.fastq.gz' -printf "%f\n" | sed 's/_R1_001.fastq.gz//'); do # remove common ending of name
 
echo "$NAME"

p1='_R1_001.fastq.gz'

echo -e '#!/bin/sh
#$ -cwd
#$ -q batchq
module add mixcr
cd /t1-data/user/lfelce/MiXCR/CD8_input/

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

cd /t1-data/user/lfelce/MiXCR/CD8_output/script/

for line in $(ls *.sh); do
sbatch $line
done

squeue -u lfelce
