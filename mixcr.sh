## below is script for mixcr. 

## run the first one as it is, it will generate individual script files for each sample

## then run second script at bottom, which will run all individual script files wherever they are saved.

cd /t1-data/user/lfelce/Dong_Pools_9-12/Data/Intensities/BaseCalls/merge/ #input directory
DIR=/t1-data/user/nassisar/SamrtSeq_Ling/ # output directory

for NAME in $(find . -name '*_R1_001.fastq.gz' -printf "%f\n" | sed 's/_R1_001.fastq.gz//'); do # remove common ending of name
 
echo "$NAME"

p1='_R1_001.fastq.gz'

echo -e '#!/bin/sh
#$ -cwd
#$ -q batchq
module add mixcr
cd /t1-data/user/lfelce/Dong_Pools_9-12/Data/Intensities/BaseCalls/merge/

mixcr align -p rna-seq -s hsa -OallowPartialAlignments=true' $NAME$p1 $DIR$NAME'.vdjca

mixcr assemblePartial' $DIR$NAME'.vdjca' $DIR$NAME'_Rescued_it1.vdjca
mixcr assemblePartial' $DIR$NAME'_Rescued_it1.vdjca' $DIR$NAME'_Rescued_it2.vdjca
mixcr assemblePartial' $DIR$NAME'_Rescued_it2.vdjca' $DIR$NAME'_Rescued.vdjca

mixcr assemble -OaddReadsCountOnClustering=true -ObadQualityThreshold=15' $DIR$NAME'_Rescued.vdjca' $DIR$NAME'.clns

mixcr exportClones' $DIR$NAME'.clns' $DIR$NAME'.txt' > $DIR'script/'$NAME'.sh'

done

#------------- run scripts on server - this should be separate script -------------------#
# cd /t1-data/user/nassisar/SamrtSeq_Ling/script/

# for line in $(ls *.sh); do
# sbatch $line
# done

# squeue -u nassisar
