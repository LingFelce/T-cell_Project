# downloading bcl files and bcl to fastq conversion
# for original file and history, see Miscellaneous repository

# basespace command line interface

# to switch user
bs auth --force

bs list runs

+---------------------------------+-----------+-----------------+----------+
|              Name               |    Id     | ExperimentName  |  Status  |
+---------------------------------+-----------+-----------------+----------+
| 200909_NB501183_0819_AHLLVLBGXG | 197260093 | Dong-Pools1-4   | Complete |
| 200910_NB501183_0820_AH7LL2BGXG | 197263145 | DongPools5-8    | Complete |
| 200910_NB502048_0489_AH7LLGBGXG | 197276079 | DongPools9-12   | Complete |
| 201104_NB501183_0841_AH2CJCBGXH | 198352159 | DONGPOOL_121020 | Complete |
+---------------------------------+-----------+-----------------+----------+

# will download all files associated with run into current folder
bs download run -i 197263145 -o ./

###### used bcl2fastq module, see below #############
# install bcl2fastq

wget ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/software/bcl2fastq/bcl2fastq2-v2-20-0-tar.zip

# use file manager mode to unzip file so that file ends .tar.gz

# as can't download rpm file, have to install from source using .tar.gz file
# (actually did manage to download rpm file, but said permission denied - presumably need admin rights)

# set up environmental variables - build and source directories must be different
export TMP=/t1-data/user/lfelce/tmp
export SOURCE=${TMP}/bcl2fastq
export BUILD=${TMP}/bcl2fastq2-build
export INSTALL_DIR=/t1-data/user/lfelce/bcl2fastq2

tar -xvzf <file name>

# followed instructions from user guide to build and install package
# https://emea.support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq2-v2-20-software-guide-15051736-03.pdf

# was getting a configuration error, realised that it was trying to use something in the conda environment tracer_env
# removed tracer_env conda environment, seems to be working fine now.

###############################################################

# convert and demultiplex BCL files - run in appropriate folder

# actually managed to get Ewan the cbrg facility manager to install bcl2fastq, so run from module instead.

# submit as job script to slurm

#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=bcl2fastq
#SBATCH --nodes=1
#SBATCH --mem=128G
#SBATCH --time=07-00:00:00
#SBATCH --output=%j_%x.out
#SBATCH --error=%j_%x.err
#SBATCH --mail-user=ling.felce@ndm.ox.ac.uk
#SBATCH --mail-type=end,fail

cd /t1-data/user/lfelce/Dong_Pools_5-8/

module load bcl2fastq

bcl2fastq --input-dir /t1-data/user/lfelce/Dong_Pools_5-8/Data/Intensities/BaseCalls/
--runfolder-dir /t1-data/user/lfelce/Dong_Pools_5-8/
--output-dir /t1-data/user/lfelce/Dong_Pools_5-8/Data/Intensities/BaseCalls/Output/
--interop-dir /t1-data/user/lfelce/Dong_Pools_5-8/Data/Intensities/BaseCalls/InterOp/
--stats-dir /t1-data/user/lfelce/Dong_Pools_5-8/Data/Intensities/BaseCalls/Stats/
--reports-dir /t1-data/user/lfelce/Dong_Pools_5-8/Data/Intensities/BaseCalls/Reports/
--sample-sheet 4x96into384demultiplex.csv 


# list of available queues and nodes
squeue

# queue info
squeue

# send job to queue
sbatch ./<jobname>.sh

# repeat with Dong_Pools_9-12
# and repeated with Dong_Pools_1-4

##########################
# submit job to SunGridEngine (BMRC cluster)

#!/bin/bash
#$ -wd /well/jknight/users/jln789/DONG231120TCR/
#$ -q short.qc
###$ -o output.log
###$ -e error.log

cd /well/jknight/users/jln789/DONG231120TCR/

module load bcl2fastq2/2.20.0-foss-2018b

bcl2fastq --input-dir /well/jknight/users/jln789/DONG231120TCR/Data/Intensities/BaseCalls/ \
--runfolder-dir /well/jknight/users/jln789/DONG231120TCR/ \
--output-dir /well/jknight/users/jln789/DONG231120TCR/Data/Intensities/BaseCalls/Output/ \
--interop-dir /well/jknight/users/jln789/DONG231120TCR/Data/Intensities/BaseCalls/InterOp/ \
--stats-dir /well/jknight/users/jln789/DONG231120TCR/Data/Intensities/BaseCalls/Stats/ \
--reports-dir /well/jknight/users/jln789/DONG231120TCR/Data/Intensities/BaseCalls/Reports/ \
--sample-sheet TCRindexing.csv

##############

# submit to queue
qsub -q short.qc ./bcl2fastq.sh
