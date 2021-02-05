# Extracting full TCR sequence from .fastq.gz files using MiXCR #

# https://mixcr.readthedocs.io/en/master/assembleContigs.html

#-----script for CCB cluster-----#

#!/bin/sh
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --mem=128G

module load mixcr

cd /t1-data/user/lfelce/MiXCR/CD8_ORF3a-28_TCR_Clones_input/

mixcr align --species hsa -p kAligner2 --report report.txt H09_S93_L001_R1_001.fastq.gz H09_S93_L001_R2_001.fastq.gz MRN345_Cl07_alignments.vdjca

# assemble default CDR3 clonotypes (note: --write-alignments is required for further contig assembly)
mixcr assemble --write-alignments --report MRN345_Cl07_report.txt MRN345_Cl07_alignments.vdjca MRN345_Cl07_clones.clna

# assemble full TCR receptors
mixcr assembleContigs --report MRN345_Cl07_report.txt MRN345_Cl07_clones.clna MRN345_Cl07_full_clones.clns

# export full TCR receptors
mixcr exportClones -c TCR -p fullImputed MRN345_Cl07_full_clones.clns MRN345_Cl07_full_clones.txt


------------------------------

# import full_clones.txt file to tidy up (remove background clones) and export as .csv
# can copy and paste target sequences into IgBLAST to double check correct TRBV/TRAV

#######################################

# https://mixcr.readthedocs.io/en/master/export.html#exporting-well-formatted-alignments-for-manual-inspection

mixcr exportAlignmentsPretty --skip 1000 --limit 10 MRN345_Cl07_alignments.vdjca test4.txt

