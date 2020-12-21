##### MiXCR code using RNA samples amplified by TCR sequencing using SMARTer Human TCR a/b Profiling Kit #####

# options
# --5-end -no-v-primers as 5'RACE with template switch oligo
# --3-end c-primers reverse primer to constant region
# --adapters no-adpaters absent from fastqc/multiqc analysis

cd /t1-data/user/lfelce/MiXCR/1131-TP-1_CD8_new_input/
DIR=/t1-data/user/lfelce/MiXCR/1131-TP-1_CD8_new_output/

for NAME in $(find . -name '*_R1_001.fastq.gz' -printf "%f\n" | sed 's/_R1_001.fastq.gz//'); do
 
echo "$NAME"

p1='_R1_001.fastq.gz'
p2='_R2_001.fastq.gz' 

echo -e '#!/bin/sh
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --mem=128G

module load mixcr

cd /t1-data/user/lfelce/MiXCR/1131-TP-1_CD8_new_output/

mixcr analyze amplicon -s hsa --starting-material rna --5-end no-v-primers --3-end c-primers --adapters no-adapters --contig-assembly --only-productive /t1-data/user/lfelce/MiXCR/1131-TP-1_CD8_new_input/'$NAME$p1 '/t1-data/user/lfelce/MiXCR/1131-TP-1_CD8_new_input/'$NAME$p2 $NAME '--receptor-type tcr' > $DIR'script/'$NAME'.sh'

done

#------------- run scripts on server
cd /t1-data/user/lfelce/MiXCR/1131-TP-1_CD8_new_output/script/

for line in $(ls 1131-TP-1_CD8_NP16_B*.sh); do
sbatch $line
done

squeue -u lfelce


