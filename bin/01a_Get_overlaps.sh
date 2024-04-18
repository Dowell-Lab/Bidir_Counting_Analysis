#!/bin/bash

#SBATCH --job-name=get_overlaps
#SBATCH --output=/scratch/Users/hoto7260/Resp_Env/Sasse2019/e_and_o/00_get_overlaps_%j.out
#SBATCH --error=/scratch/Users/hoto7260/Resp_Env/Sasse2019/e_and_o/00_get_overlaps_%j.out
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=5Gs
#SBATCH --partition short
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=hoto7260@colorado.edu

module load bedtools
##########################
# EDIT THE FOLLOWING
##########################
### Naming
WD=
OUT_DIR=${WD}/overlaps
PREFIX="Sasse_2019_nascent"
DATE=04-4-24
### Files
# muMerge Bid file for counting (I have used 600bp window)
COUNT_WIN=/Users/hoto7260/projects/Resp_Env/Sasse2019/data/beds/${PREFIX}_${DATE}_MUMERGE_600bpwin_tfit,dreg.sorted.bed
# MuMerge Bid file for TSS identification only (I have used 50bp window)
TSS_WIN=/Users/hoto7260/projects/Resp_Env/Sasse2019/data/beds/${PREFIX}_${DATE}_MUMERGE_50bpwin_tfit,dreg.sorted.bed
# The FULL transcripts with putatives (WARNING: if names are too long, bedtools breaks without warning)
TRANSCRIPTS=/scratch/Shares/dowell/hoto7260/Bidir_Count/hg38_refseq_diff53prime_with_putatives_fixnames.bed
# The gene file you will use for counting (make sure truncated)
TRUNC=/scratch/Shares/dowell/genomes/hg38/ncbi/hg38_refseq_diff53prime_5ptrunc_counting.bed
# TSS 1kb region (500bp ± TSS)
TSS_BED=/scratch/Shares/dowell/hoto7260/Bidir_Count/hg38_TSS1kb_refseq_diff53prime_with_putatives.bed

##########################
# NAMING THE OUTPUT FILES
##########################
mkdir -p ${OUT_DIR}
TSS_OUT=${OUT_DIR}/overlaps_hg38_TSS1kb_withput_${PREFIX}_MUMERGE_${DATE}.bed
COUNT_OUT=${OUT_DIR}/overlaps_hg38_withput_${PREFIX}_MUMERGE_${DATE}.bed
COUNT_TRUNC_OUT=${OUT_DIR}/overlaps_hg38_trunc_${PREFIX}_MUMERGE_${DATE}.bed
CLOSE_OUT=${OUT_DIR}/closest_hg38_withput_${PREFIX}_MUMERGE_${DATE}.bed

##########################
# GETTING THE OVERLAPS
##########################
# Identifying TSS Bids
bedtools intersect -wo -a ${TSS_WIN} -b ${TSS_BED} > ${TSS_OUT}
wc -l ${TSS_OUT}

## Bids that overlap with transcripts (for counting) -- must have gene first
bedtools intersect -wo -a ${TRANSCRIPTS} -b ${COUNT_WIN} > ${COUNT_OUT}
wc -l ${COUNT_OUT}

## Bids that overlap with transcripts (for counting) -- must have gene first
bedtools intersect -wo -a ${TRUNC} -b ${COUNT_WIN} > ${COUNT_TRUNC_OUT}
wc -l ${COUNT_TRUNC_OUT}

## Bids within 10kb of the transcripts (do 1000 to ensure all are included)
bedtools closest -D ref -k 1000 -a ${TRANSCRIPTS} -b ${COUNT_WIN} > ${CLOSE_OUT}
wc -l ${CLOSE_OUT}