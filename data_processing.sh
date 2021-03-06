# merging lanes for .fastq.gz files

# fine unique identifier for each file
ls -1 *_R1*.gz | awk -F '_' '{print $1"_"$2"_"$3}' | sort | uniq > ID

# merge lanes L001, L002, L003, L004 - only need to do once as single end reads, no R2
for i in `cat ./ID`; do cat $i\_L001_R1_001.fastq.gz  $i\_L002_R1_001.fastq.gz  $i\_L003_R1_001.fastq.gz  $i\_L004_R1_001.fastq.gz > merge/$i\_R1_001.fastq.gz; done

# put in batch script and send to queue (can do ID bit separately)


# check md5sum
md5sum *fastq.gz > md5sum.txt
md5sum -c md5sum.txt

# rename files
# can use command line for i in *; do mv "$i" 1_"$i"; done but need to specify i first
# otherwise can go to File Manager, select multiple files and rename - can also preview new name before making changes

# soft links to wherever going to use .fastq.gz files
find /t1-data/user/lfelce/Dong_Pools_5-8/Data/Intensities/BaseCalls/merge -name "022_*.fastq.gz" | xargs -I v_f ln -s v_f
find /t1-data/user/lfelce/Dong_Pools_5-8/Data/Intensities/BaseCalls/merge -name "1062_*.fastq.gz" | xargs -I v_f ln -s v_f
find /t1-data/user/lfelce/Dong_Pools_9-12/Data/Intensities/BaseCalls/merge -name "025_*.fastq.gz" | xargs -I v_f ln -s v_f
find /t1-data/user/lfelce/Dong_Pools_9-12/Data/Intensities/BaseCalls/merge -name "1504_*.fastq.gz" | xargs -I v_f ln -s v_f
find /t1-data/user/lfelce/Dong_Pools_9-12/Data/Intensities/BaseCalls/merge -name "1493_*.fastq.gz" | xargs -I v_f ln -s v_f
find /t1-data/user/lfelce/Dong_Pools_9-12/Data/Intensities/BaseCalls/merge -name "1525_*.fastq.gz" | xargs -I v_f ln -s v_f

# CD8 samples
find /t1-data/user/lfelce/Dong_Pools_5-8/Data/Intensities/BaseCalls/merge -name "1131-*.fastq.gz" | xargs -I v_f ln -s v_f
find /t1-data/user/lfelce/Dong_Pools_5-8/Data/Intensities/BaseCalls/merge -name "1201-*.fastq.gz" | xargs -I v_f ln -s v_f

find /t1-data/user/lfelce/Dong_Pools_1-4/Data/Intensities/BaseCalls/merge -name "005_*.fastq.gz" | xargs -I v_f ln -s v_f
find /t1-data/user/lfelce/Dong_Pools_1-4/Data/Intensities/BaseCalls/merge -name "1153_*.fastq.gz" | xargs -I v_f ln -s v_f
find /t1-data/user/lfelce/Dong_Pools_1-4/Data/Intensities/BaseCalls/merge -name "1131-*.fastq.gz" | xargs -I v_f ln -s v_f

