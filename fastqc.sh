# code for running fastqc in queue
# make sure create output folder first, easier to do in same folder where .fastq files are stored

#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=fastqc
#SBATCH --nodes=1
#SBATCH --mem=128G
#SBATCH --time=07-00:00:00
#SBATCH --output=%j_%x.out
#SBATCH --error=%j_%x.err

cd /t1-data/user/lfelce/MiXCR/CD4_input

module load fastqc

for i in /t1-data/user/lfelce/MiXCR/CD4_input/*.fastq.gz
do
  fastqc -q -t 12 --nogroup ${i} --outdir fastqc
done
