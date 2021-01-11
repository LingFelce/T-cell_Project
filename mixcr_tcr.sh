##### MiXCR code using RNA samples amplified by TCR sequencing using SMARTer Human TCR a/b Profiling Kit #####

# options
# --5-end -no-v-primers as 5'RACE with template switch oligo
# --3-end c-primers reverse primer to constant region
# --adapters no-adpaters absent from fastqc/multiqc analysis

cd /t1-data/user/lfelce/MiXCR/CD8_ORF3a-28_TCR_Clones_input/
DIR=/t1-data/user/lfelce/MiXCR/CD8_ORF3a-28_TCR_Clones_output/

for NAME in $(find . -name '*_R1_001.fastq.gz' -printf "%f\n" | sed 's/_R1_001.fastq.gz//'); do
 
echo "$NAME"

p1='_R1_001.fastq.gz'
p2='_R2_001.fastq.gz' 

echo -e '#!/bin/sh
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --mem=128G

module load mixcr

cd /t1-data/user/lfelce/MiXCR/CD8_ORF3a-28_TCR_Clones_output/

mixcr analyze amplicon -s hsa \
--starting-material rna \
--5-end no-v-primers \
--3-end c-primers \
--adapters no-adapters \
--contig-assembly \
--only-productive /t1-data/user/lfelce/MiXCR/CD8_ORF3a-28_TCR_Clones_input/'$NAME$p1 '/t1-data/user/lfelce/MiXCR/CD8_ORF3a-28_TCR_Clones_input/'$NAME$p2 \
$NAME '--receptor-type tcr' > $DIR'script/'$NAME'.sh'

done


#------------- run scripts on server
cd /t1-data/user/lfelce/MiXCR/CD8_ORF3a-28_TCR_Clones_output/script/

for line in $(ls B*.sh); do
sbatch $line
done

squeue -u lfelce

#####################################
# run mixcr on BMRC cluster
# have to specify separate receptor type tra and trb in two separate scripts (run in separate folders)
# otherwise TRB.txt file was same as TRA.txt file!

cd /well/jknight/users/jln789/TCR071220_Plate-2/Data/Intensities/BaseCalls/Output/
DIR=/well/jknight/users/jln789/mixcr/TCR071220_Plate-2/

for NAME in $(find . -name '*_R1_001.fastq.gz' -printf "%f\n" | sed 's/_R1_001.fastq.gz//'); do
 
echo "$NAME"

p1='_R1_001.fastq.gz'
p2='_R2_001.fastq.gz' 

echo -e '#!/bin/bash
#$ -wd /well/jknight/users/jln789/mixcr/TCR071220_Plate-2/
#$ -q short.qc

module load MiXCR/3.0.3-Java-1.8

cd /well/jknight/users/jln789/mixcr/TCR071220_Plate-2/

mixcr analyze amplicon -s hsa \
--starting-material rna \
--5-end no-v-primers \
--3-end c-primers \
--adapters no-adapters \
--contig-assembly \
--only-productive /well/jknight/users/jln789/TCR071220_Plate-2/Data/Intensities/BaseCalls/Output/'$NAME$p1 '/well/jknight/users/jln789/TCR071220_Plate-2/Data/Intensities/BaseCalls/Output/'$NAME$p2 \
$NAME '--receptor-type trb' > $DIR'script/'$NAME'.sh'

done

######

# submit individual scripts to queue

cd /well/jknight/users/jln789/mixcr/TCR071220_Plate-2/script

for line in $(ls *.sh); do
qsub -q short.qc $line
done

qstat
