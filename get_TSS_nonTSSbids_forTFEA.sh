#!/bin/bash

#SBATCH --job-name=Bidir_Counting
#SBATCH --output=/scratch/Users/path/e_and_o/comb_get_counts_%j.out
#SBATCH --error=/scratch/Users/path/e_and_o/comb_get_counts_%j.err
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=5G
#SBATCH --partition short
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=hoto7260@colorado.edu

##########################
# EDIT THE FOLLOWING
##########################

#### Loading Modules/Environments
# if you have a conda or other virual environment, activate here
source activate /Users/hoto7260/miniconda3/envs/Rjupyter (in my case, this environment has  R with data.table & stringR AND python 3.9 with pandas and numpy)

module load bedtools
module load samtools/1.8
module load subread/1.6.2
module load R/4.4.0

### Paths
SRC=/Users/hoto7260/src/Bidir_Counting_Analysis/
WD=/scratch/Users/hoto7260/Bench_DE/Elkon2015myc
### Naming
PREFIX="Elkon2015myc"
DATE=04_4_24

## IN Directories (BAMS in this case)
bams=/scratch/Users/hoto7260/Bench_DE/raw_data/bams/Elkon2015myc/PRO/counting/

### Files
File with consensus regions (e.g. output from Mumerge)
CONS_FILE=

### Parameters
# window for counting (mu-window & mu + window)
COUNT_WIN=500
# window for TSS overlaps (mu-window & mu + window)
TSS_WIN=25
# minimum rowSUM counts needed for a gene to be considered "transcribed" and significantly altering bid counts (I used 40 for 4 samples)
COUNT_LIMIT_GENES=60
# the number required for a bidirectional to be considered impeding on gene transcription
# usually do 25 for 4 samples
COUNT_LIMIT_BIDS=30
# do you want to get TFEA combined file outputs? YES or NO
TFEA="YES"
# do you want to have unique names be created for each bidirectional? YES or NO
MAKE_NAMES="YES"

echo "##########################################################################"
echo "########################            VARIABLES         ####################"
echo "##########################################################################"
echo "Output Directory" ${WD}
echo "Prefix" ${PREFIX}
echo "Date" ${DATE}
echo "Consensus file for bids" ${COUNT_WIN}
echo "Window used for counting (*2 for total region width)" ${COUNT_WIN}
echo "Window used for TSS overlaps (*2 for total region width)" ${TSS_WIN}
echo "min rowSum counts for gene to be a convolution problem" ${COUNT_LIMIT_GENES}
echo "min rowSum counts for a bid to be a convolution problem" ${COUNT_LIMIT_BIDS}
echo ""


echo "#####################################################################################"
echo "########################            1a. GET THE OVERLAPS         ####################"
echo "#####################################################################################"

# The FULL transcripts with putatives (WARNING: if names are too long, bedtools breaks without warning)
TRANSCRIPTS=${SRC}/assets/hg38_refseq_diff53prime_with_putatives_fixnames.bed
# The gene file you will use for counting (make sure truncated)
TRUNC=${SRC}/assets/hg38_refseq_diff53prime_5ptrunc_counting.bed
# TSS 1kb region (500bp ± TSS)
TSS_BED=${SRC}/assets/hg38_TSS1kb_refseq_diff53prime_with_putatives.bed

##########################
# GET THE INPUT BED FILES
# mumerge_filename, count_window, tss_window, prefix_date
Rscript ${SRC}/bin/get_window_cons_files.r ${CONS_FILE} ${COUNT_WIN} ${TSS_WIN} ${WD}/regions/${PREFIX}_${DATE} ${MAKE_NAMES}    
COUNT_WIN_FILE=${WD}/regions/${PREFIX}_${DATE}_MUMERGE_${COUNT_WIN}win_count.sorted.bed
TSS_WIN_FILE=${WD}/regions/${PREFIX}_${DATE}_MUMERGE_${COUNT_WIN}win_TSS.sorted.bed
##########################
# NAMING THE OUTPUT FILES
OUT_DIR=${WD}/overlaps
mkdir -p ${OUT_DIR}
TSS_OUT=${OUT_DIR}/overlaps_hg38_TSS1kb_withput_${PREFIX}_MUMERGE_${DATE}.bed
COUNT_OUT=${OUT_DIR}/overlaps_hg38_withput_${PREFIX}_MUMERGE_${DATE}.bed
COUNT_TRUNC_OUT=${OUT_DIR}/overlaps_hg38_trunc_${PREFIX}_MUMERGE_${DATE}.bed
CLOSE_OUT=${OUT_DIR}/closest_hg38_withput_${PREFIX}_MUMERGE_${DATE}.bed

# ##########################
# # GETTING THE OVERLAPS
# first make sure I'm only using the first four columns
INT_COUNT_WIN=${WD}/int_count_win.bed
cut -f1-4 ${COUNT_WIN_FILE} > ${INT_COUNT_WIN}
COUNT_WIN_FILE=${INT_COUNT_WIN}
INT_TSS_WIN=${WD}/int_tss_win.bed
cut -f1-4 ${TSS_WIN_FILE} > ${INT_TSS_WIN}
TSS_WIN_FILE=${INT_TSS_WIN}

# Identifying TSS Bids
bedtools intersect -wo -a ${TSS_WIN_FILE} -b ${TSS_BED} > ${TSS_OUT}
wc -l ${TSS_OUT}
## Bids that overlap with transcripts (for counting) -- must have gene first
bedtools intersect -wo -a ${TRANSCRIPTS} -b ${COUNT_WIN_FILE} > ${COUNT_OUT}
wc -l ${COUNT_OUT}
## Bids that overlap with transcripts (for counting) -- must have gene first
bedtools intersect -wo -a ${TRUNC} -b ${COUNT_WIN_FILE} > ${COUNT_TRUNC_OUT}
wc -l ${COUNT_TRUNC_OUT}
## Bids within 10kb of the transcripts (do 1000 to ensure all are included)
bedtools closest -D ref -k 100 -a ${TRANSCRIPTS} -b ${COUNT_WIN_FILE} > ${CLOSE_OUT}
wc -l ${CLOSE_OUT}

echo ""
echo "#####################################################################################"
echo "########################          1b. GET PREL GENE COUNTS       ####################"
echo "#####################################################################################"
counts=${WD}/counts
fixed_counts=${WD}/fixed_counts
mkdir -p ${counts}
mkdir -p ${fixed_counts}

genes=${SRC}/assets/hg38_refseq_diff53prime_with_putatives_5ptrunc.saf

########################################## 
# Count reads over gene coordinates     ##
##########################################
cd ${bams}
echo "Submitted counts with Rsubread for original genes by strand......"
# Use this many threads (-T)
# Count multi-overlapping read
# Count reads in a strand specific manner (-s 2)
# Count by the gene feature (-t 'gene')
featureCounts \
    -T 32 \
    -O \
    -s 1 \
    -a ${genes} \
    -F 'SAF' \
    -o ${counts}/${PREFIX}_str_put_genes.txt \
    ./*.sorted.bam
    
echo ""
echo "###################################################################################"
echo "########################          2. GET TSS & Filt Bids       ####################"
echo "###################################################################################"
# YES refers to TFEA output being made
Rscript ${SRC}/bin/get_TSS_and_filt_bids.r ${WD} ${PREFIX} ${DATE} ${COUNT_WIN_FILE} ${COUNT_LIMIT_GENES} ${TFEA}
wc -l ${WD}/regions/*

echo ""
echo "###################################################################################"
echo "########################         3. GET Fixed & Bid Counts     ####################"
echo "###################################################################################"
bidirs=${WD}/regions/${PREFIX}_MUMERGE_tfit,dreg_${DATE}_filt.saf

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
   -o ${counts}/${PREFIX}_uns_bidirs.txt \
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
   -o ${counts}/${PREFIX}_pos_bidirs.txt \
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
   -o ${counts}/${PREFIX}_neg_bidirs.txt \
   ./*.sorted.bam

###################
# REMOVE INTERMEDIATE FILES #
###################
rm ${TSS_WIN_FILE}
rm ${COUNT_WIN_FILE}


echo "DONE!"
date
