#!/bin/bash

#SBATCH --job-name=featureCounts_bids_SE
#SBATCH --output=/scratch/Users/hoto7260/Resp_Env/Sasse2019/e_and_o/featurecounts_bids_%j.out
#SBATCH --error=/scratch/Users/hoto7260/Resp_Env/Sasse2019/e_and_o/featurecounts_bids_%j.err
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=5G
#SBATCH --partition short
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hoto7260@colorado.edu

##load modules
module load samtools/1.8
module load subread/1.6.2

##########################
# EDIT THE FOLLOWING
##########################
## Naming
prefix="Sasse2019_nascent"
date=04_04_24
## IN/OUT Directories
# wd is where the bams and counts will be saved
wd=/scratch/Users/hoto7260/Resp_Env/Sasse2019/out
bams=${wd}/bams
counts=${wd}/counts
fixed_counts=${wd}/fixed_counts 

bidirs=${wd}/regions/${prefix}_MUMERGE_tfit,dreg_${date}_filt.saf

cd ${bams}

echo "Submitted counts with Rsubread for bidsirectionals unstranded......"
# Use this many threads
# Count multi-overlapping read (regardless of percentage of overlap)
# Count reads NOT in a strand specific manner
featureCounts \
   -T 32 \
   -O \
   -s 0 \
   -a ${bidirs} \
   -F 'SAF' \
   -o ${counts}/${prefix}_uns_bidirs.txt \
   ./*.sorted.bam

echo "Submitted counts with Rsubread for bidirectionals stranded (+)......"
# Use this many threads
# Count multi-overlapping read
# Count reads NOT in a strand specific manner
featureCounts \
   -T 32 \
   -O \
   -s 1 \
   -a ${bidirs} \
   -F 'SAF' \
   -o ${counts}/${prefix}_pos_bidirs.txt \
   ./*.sorted.bam

echo "Submitted counts with Rsubread for bidirectionals stranded (-)......"
# Use this many threads
# Count multi-overlapping read
# Count reads NOT in a strand specific manner
featureCounts \
   -T 32 \
   -O \
   -s 2 \
   -a ${bidirs} \
   -F 'SAF' \
   -o ${counts}/${prefix}_neg_bidirs.txt \
   ./*.sorted.bam
