# Processing and preparing data for single cell RNA-Seq analysis - from .fastq.gz to counts matrix

# Check on a few cells first - maybe 1 from each patient? CD8 or CD4?

# check if have appropriate modules installed

# .fastq.gz quality control

mkdir fastqc_results
fastqc -o fastqc_results <input file>

# trimming - necessary?

# find reference transcriptome to for STAR - reference genome sequences (FASTA) and annotations (GTF) - produce genome index

mkdir indices
mkdir indices/STAR
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir indices/STAR --genomeFastaFiles <fastafile> --sjdbGTFfile <GTFfile>

# mapping using STAR - output directly as BAM
# read files are compressed --readFilesCommand zcat or --readFilesCommand gunzip -c
# map multiple samples at same time - use commas to separate single-end reads with no space 
# output as BAM so don't need to convert SAM files to BAM - check if need unsorted or sorted by coordinates (as input for featureCounts)

mkdir results
mkdir results/STAR

STAR --runThreadN 4 --genomeDir indices/STAR --readFilesIn <fastqfile1>,<fastqfile2>,<fastqfile3> --outFileNamePrefix results/STAR/prefix

# mapping quality control - multiqc

# use featureCounts to quantify reads (check what multimapping is again!)

# include multimapping
<featureCounts_path>/featureCounts -O -M -Q 30 -p -a genome.gtf -o outputfile input.bam

# exclude multimapping
<featureCounts_path>/featureCounts -Q 30 -p -a genome.gtf -o outputfile input.bam

featureCounts -t exon -g gene_id -a <GTF> -o counts.txt <bam1> <bam2> <bam3>

# once have counts matrix, can analyse in R





