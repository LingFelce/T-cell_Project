##### Cell Ranger pipeline for analysing 10x data #####

## bcl conversion to fastq ##

cellranger mkfastq --help

Run Illumina demultiplexer on sample sheets that contain 10x-specific sample 
index sets, and generate 10x-specific quality metrics after the demultiplex.  
Any bcl2fastq argument works, except a few that are set by the pipeline 
to ensure proper trimming and sample indexing. The FASTQ output generated 
is the same as when running bcl2fastq directly.

cellranger mkfastq --id=tutorial_walk_through \
--run=/mnt/home/user.name/yard/run_cellrnager_mkfastq/cellranger-tiny-bcl-1.2.0 \
--csv=/mnt/home/user.name/yard/run_cellrnager_mkfastq/cellranger-tiny-bcl-simple-1.2.0.csv
