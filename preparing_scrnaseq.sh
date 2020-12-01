# Processing and preparing data for single cell RNA-Seq analysis - from .fastq.gz to counts matrix

# Check on a few cells first - maybe 1 from each patient? CD8 or CD4?

# check if have appropriate modules installed

# .fastq.gz quality control - see separate script as send to queue, module available

mkdir fastqc_results
fastqc -o fastqc_results <input file>

# multiqc installed on conda so have to run locally, but doesn't take very long

# trimming - necessary? if so, use trim_galore
# no adapter content, so don't need to trim

# find reference transcriptome to for STAR - reference genome sequences (FASTA) and annotations (GTF) - produce genome index

# don't need to make genome index, already have on cbrg cluster /databank/indices/star/hg19

# mapping using STAR - output directly as BAM
# read files are compressed --readFilesCommand zcat or --readFilesCommand gunzip -c
# map multiple samples at same time - use commas to separate single-end reads with no space 
# output as BAM so don't need to convert SAM files to BAM - check if need unsorted or sorted by coordinates (as input for featureCounts)

# softlink .fastq.gz files from MiXCR input folders to STAR folder
find /t1-data/user/lfelce/MiXCR/CD4_input -name "*.fastq.gz" | xargs -I v_f ln -s v_f

#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=star_map
#SBATCH --nodes=1
#SBATCH --mem=128G
#SBATCH --time=07-00:00:00
#SBATCH --output=%j_%x.out
#SBATCH --error=%j_%x.err

module load STAR/2.6.1d
cd /t1-data/user/lfelce/STAR/

for i in ./*_R1_001.fastq.gz

do

STAR --runThreadN 10 --readFilesCommand gunzip -c --outSAMtype BAM SortedByCoordinate --genomeDir /databank/indices/star/hg19 --readFilesIn $i --outFileNamePrefix ./results/${i%_R1_001.fastq.gz}_

done


# mapping quality control - multiqc

# use featureCounts to quantify reads 
# featureCounts is part of subreads package
# exclude multimapping (10x Cell Ranger by default exludes mutimapping, so must be ok!)

<featureCounts_path>/featureCounts -Q 30 -p -a genome.gtf -o outputfile input.bam

featureCounts -t exon -g gene_id -a <GTF> -o counts.txt <bam1> <bam2> <bam3> #-t and -g are default, -Q minimum mapping score 0 is default



# once have counts matrix, can analyse in R - use R on cluster?

# modules available
alicevision/2.2.0              bowtie2/2.3.4.3                 fastqc/0.11.9            java/7u80                         python-cbrg/202008            seqtk/20201102          
allelecounter-castel/20181017  bowtie2/2.4.1                   fiji/20190409            java/8u131                        python-cbrg/202009            snpeff/5.0              
annovar/20180416               bpipe/0.9.9.9(default)          flair/1.5                java/8u171                        python-cbrg/202010            snptest/2.5.2           
annovar/20191024               brie/0.2.0                      gatk/4.0.5.2             java/8u191                        python-cbrg/202011            snver/0.5.3             
aspera/3.9.8                   bwa/0.7.17                      gatk/4.0.11.0            julia/1.4.1                       python-cbrg/current(default)  sratoolkit/2.10.8       
aspera/3.9.9                   catt/1.8                        gcc/8.2.0                libBigWig/0.4.4                   qctool/2.0.7                  stampy/1.0.32           
atacqc/20170728                cbrg/centos6                    gcc/9.3.0                manta/1.6.0                       qtltools/1.2                  STAR/2.6.0c             
backspin/1.0(default)          cbrg/centos7(default)           gmp/6.2.0                mathematica/11.3.0                R-base/3.6.3                  STAR/2.6.1d             
bam-readcount/20201113         ccs/4.2.0                       gsea/3.0                 matlab/r2018b                     R-base/4.0.1(default)         STAR/2.7.3a(default)    
bamscale/0.0.5                 cellpose/20200526               gsl/2.6                  matlab/r2019b                     R-cbrg/202006                 stress/1.0.4            
bamscale/1.0                   cellprofiler/3.1.9              gtool/0.7.5              matlab/r2020a                     R-cbrg/202007                 subread/2.0.0           
bamtools/2.5.1                 cellranger/2.1.1                guppy-gpu/3.3.3          meme/5.1.1(default)               R-cbrg/202008                 tensorflow1-cpu/202010  
bamutil/1.0.14                 cellranger/3.0.1                guppy/3.3.3              meshroom/2019.2.0                 R-cbrg/202009                 tensorflow1-gpu/202010  
bcftools/1.5                   cellranger/3.1.0                hdf5/1.10.7              minimap2/2.17                     R-cbrg/202010                 tensorflow2-cpu/202010  
bcftools/1.9                   cellranger/4.0.0                hdfview/3.1.1            mixcr/3.0.13                      R-cbrg/202011                 tensorflow2-gpu/202010  
bcftools/1.10.2                circexplorer/1.1.10             hisat2/20201102          moods/1.9.4.1                     R-cbrg/current(default)       tophat/2.1.1            
bcl2fastq/2.20.0.422           clairvoyante/1.02               htslib/1.5               mpileup2readcounts/1.0            rna-star/2.6.0c               trim_galore/0.6.5       
bedtools/2.26.0                cmake/3.17.2                    htslib/1.9               mpileup2readcounts/1.0a(default)  rna-star/2.6.1d               ucsctools/369           
bedtools/2.27.1                conifer/0.2.2(default)          htslib/1.10              octopus/20200902                  rna-star/2.7.3a(default)      ucsctools/385           
bedtools/2.29.2                cuda/8.0                        htslib/1.10.2(default)   omero/5.2.4                       rstudio/1.0.143               unet/20190409           
bedtools2/2.26.0               cuda/9.0                        igblast/1.16.0           perl/5.32.0                       rstudio/1.1.383               varscan/2.4.2           
bedtools2/2.27.1               cuda/10.0                       igor/1.4.0               picard-tools/2.18.17              rstudio/1.2.5042              vcflib/20201113         
bedtools2/2.29.2               cuda/10.1(default)              igv/2.4.11               pigz/2.4                          samblaster/0.1.24             wxpython/4.0.7post2     
BETA/1.0.7(default)            cutadapt/2.10                   ilastik/1.3.3post3       pycharm/2020.1.3                  samplot/1.0.1                 zifa/20180303(default)  
beta/1.0.7(default)            cytoscape/3.6.1                 imacyte/20200929         python-base/2.7.18                samtools/1.5                  
boost/1.74.0                   discovar-denovo/52488(default)  imctools/1.0.8(default)  python-base/3.6.10                samtools/1.9                  
bowtie/1.2.1.1                 discovar/52488(default)         imctools/2.1.0           python-base/3.8.3(default)        samtools/1.10(default)        
bowtie/1.2.3(default)          expedition/20201112             impute2/2.3.2            python-cbrg/202006                scurgen/20180711              
bowtie2/2.3.2                  fastqc/0.11.5                   jags/4.3.0               python-cbrg/202007                seqkit/0.12.1


