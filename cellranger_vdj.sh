# Optimisation of TCR sequencing from Peng's data - Dong231120 
# Patient 1131-TP-1 single cell Smart Seq RNA-Seq and targeted TCR single cell sequencing (separate runs)
# Targeted TCR sequencing using SMARTer Human TCR a/b Profiling Kit
# Tried MiXCR mixcr analyze amplicon but gave too many alpha and beta chains per cell. 
# Also tried MiXCR mixcr analyze shotgun on SmartSeq2 fastq files but didn't have enough MiXCR output for single cells, only bulk

# example from 10X website

cellranger vdj --id=sample345 \ # unique run ID strong
                 --reference=/opt/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0 \ # path to CellRanger compatible reference
                 --fastqs=/home/jdoe/runs/HAWT7ADXX/outs/fastq_path \ # path of FASTQ folder - can take multiple comma-separated paths but will treat all reads from library as one sample
                 --sample=mysample \ # sample name as specified in sample sheet supplied to mkfastq. Can take multiple comma-separated values, if sample name and fastq file prefix not identical
                 --localcores=8 \ # restricts cellranger to use specified number of cores
                 --localmem=64 # restricts cell ranger to specified amount of memory



# From Charlotte's pipeline
for sample in ${samples[@]}; do

cmd="cellranger vdj --id{sample} --sample=${sample} --fastqs=${fastqs} --reference=${ref} --jobmode=${mode}

done

samples="TuPo TuNg NoPo NoNg"
fastqs="outs/fastq_path/"
transcriptome="refdata-cellranger-vdj-GrCh38-alts-ensembl-2.0.0"
ref="/databank/10x-rangers/${transcriptome}"



cellranger vdj --id=
