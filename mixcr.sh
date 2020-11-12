# batch script for mixcr

#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=mixcr_align
#SBATCH --nodes=1
#SBATCH --mem=128G
#SBATCH --time=07-00:00:00
#SBATCH --output=%j_%x.out
#SBATCH --error=%j_%x.err
#SBATCH --mail-user=ling.felce@ndm.ox.ac.uk
#SBATCH --mail-type=end

cd /t1-data/user/lfelce/MiXCR/

module load mixcr

for index in /t1-data/user/lfelce/MiXCR/*.fastq.gz
do
  mixcr align -p rna-seq -s hsa -r ${index}_alignmentReport.txt -OallowPartialAlignments=true ${index} ${index}_alignments.vdjca
done

# seems to be working

# alignment for each file taking a long time and can't assign more than 1 node for some reason, so made 6 different batch scripts to separate out samples
# 022_*.fastq.gz, 1062_*.fastq.gz, 025_*.fastq.gz, 1504_*.fastq.gz, 1493_*.fastq.gz and 1525_*.fastq.gz
# sent 6 batch scripts to queue separately, equivalent of running 6 nodes

# follow https://mixcr.readthedocs.io/en/master/rnaseq.html and/or https://mixcr.readthedocs.io/en/master/quickstart.html#analysis-of-rna-seq-data

# assemble partial reads - need to do twice

#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=mixcr_partial
#SBATCH --nodes=1
#SBATCH --mem=128G
#SBATCH --time=07-00:00:00
#SBATCH --output=%j_%x.out
#SBATCH --error=%j_%x.err
#SBATCH --mail-user=ling.felce@ndm.ox.ac.uk
#SBATCH --mail-type=end

cd /t1-data/user/lfelce/MiXCR/

module load mixcr

for index in /t1-data/user/lfelce/MiXCR/*.vdjca
do
  mixcr assemblePartial ${index} ${index}_rescued.vdjca
done

#########################

for index in /t1-data/user/lfelce/MiXCR/*_rescued.vdjca
do
  mixcr assemblePartial ${index} ${index}_2.vdjca
done

