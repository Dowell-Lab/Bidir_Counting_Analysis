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

source activate /Users/hoto7260/miniconda3/envs/Rjupyter (must have R with data.table & stringR AND python 3.9 with pandas and numpy)

module load bedtools
module load samtools/1.8
module load subread/1.6.2

##########################
# EDIT THE FOLLOWING
##########################
SRC=/Users/hoto7260/src/Bidir_Counting_Analysis/bin
### Naming
WD=/scratch/Users/hoto7260/Bench_DE/Elkon2015myc
OUT_DIR=${WD}/overlaps
PREFIX="Elkon2015myc"
DATE=04_4_24

## IN Directories (BAMS in this case)
bams=/scratch/Users/hoto7260/Bench_DE/raw_data/bams/Elkon2015myc/PRO/counting/

### Files
# muMerge Bid file for counting (I have used 600bp window)
COUNT_WIN=${WD}/regions/${PREFIX}_${DATE}_MUMERGE_600bpwin_tfit,dreg.sorted.bed
# MuMerge Bid file for TSS identification only (I have used 50bp window)
TSS_WIN=${WD}/regions/${PREFIX}_${DATE}_MUMERGE_50bpwin_tfit,dreg.sorted.bed

# minimum rowSUM counts needed for a gene to be considered "transcribed" and significantly altering bid counts (I used 40 for 4 samples)
COUNT_LIMIT_GENES=60
# the number required for a bidirectional to be considered impeding on gene transcription
# usually do 25 for 4 samples
COUNT_LIMIT_BIDS=30

echo "##########################################################################"
echo "########################            VARIABLES         ####################"
echo "##########################################################################"
echo "Output Directory" ${WD}
echo "Prefix" ${PREFIX}
echo "Date" ${DATE}
echo "Bid file for counting" ${COUNT_WIN}
echo "Bid file for identifying TSS" ${TSS_WIN}
echo "min rowSum counts for gene to be a convolution problem" ${COUNT_LIMIT_GENES}
echo "min rowSum counts for a bid to be a convolution problem" ${COUNT_LIMIT_BIDS}
echo ""


echo "#####################################################################################"
echo "########################            1a. GET THE OVERLAPS         ####################"
echo "#####################################################################################"

# The FULL transcripts with putatives (WARNING: if names are too long, bedtools breaks without warning)
TRANSCRIPTS=/scratch/Shares/dowell/hoto7260/Bidir_Count/hg38_refseq_diff53prime_with_putatives_fixnames.bed
# The gene file you will use for counting (make sure truncated)
TRUNC=/scratch/Shares/dowell/genomes/hg38/ncbi/hg38_refseq_diff53prime_5ptrunc_counting.bed
# TSS 1kb region (500bp ± TSS)
TSS_BED=/scratch/Shares/dowell/hoto7260/Bidir_Count/hg38_TSS1kb_refseq_diff53prime_with_putatives.bed

##########################
# NAMING THE OUTPUT FILES

mkdir -p ${OUT_DIR}
TSS_OUT=${OUT_DIR}/overlaps_hg38_TSS1kb_withput_${PREFIX}_MUMERGE_${DATE}.bed
COUNT_OUT=${OUT_DIR}/overlaps_hg38_withput_${PREFIX}_MUMERGE_${DATE}.bed
COUNT_TRUNC_OUT=${OUT_DIR}/overlaps_hg38_trunc_${PREFIX}_MUMERGE_${DATE}.bed
CLOSE_OUT=${OUT_DIR}/closest_hg38_withput_${PREFIX}_MUMERGE_${DATE}.bed

# ##########################
# # GETTING THE OVERLAPS

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
bedtools closest -D ref -k 100 -a ${TRANSCRIPTS} -b ${COUNT_WIN} > ${CLOSE_OUT}
wc -l ${CLOSE_OUT}

echo ""
echo "#####################################################################################"
echo "########################          1b. GET PREL GENE COUNTS       ####################"
echo "#####################################################################################"
counts=${WD}/counts
fixed_counts=${WD}/fixed_counts
mkdir -p ${counts}
mkdir -p ${fixed_counts}

genes=/scratch/Shares/dowell/hoto7260/Bidir_Count/hg38_refseq_diff53prime_with_putatives_5ptrunc.saf

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

Rscript ${SRC}/get_TSS_and_filt_bids.r ${WD} ${PREFIX} ${DATE} ${COUNT_WIN} ${COUNT_LIMIT_GENES}
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

#########################
## Now Fix the Counts ##
#########################
Rscript ${SRC}/fix_bid_counts.r ${WD} ${PREFIX} ${DATE} ${COUNT_LIMIT_BIDS} ${COUNT_LIMIT_GENES}
wc -l ${fixed_counts}/*

echo ""
echo "###################################################################################"
echo "########################         3. GET Fixed & Gene Counts     ####################"
echo "###################################################################################"

###################################
# GET THE FIXED REGIONS FOR GENES #
###################################
USE_BID=${WD}/regions/${PREFIX}_nontssbid_uniqueid_above${COUNT_LIMIT_BIDS}_${DATE}.bed
NEW_REG=${WD}/regions/${PREFIX}_new_regions_len_posneg_gene_filt_divconv_diff53prim_5ptrunc_${DATE}.gtf
COUNT_STATS=${WD}/regions/${PREFIX}_stats_len_posneg_gene_filt_divconv_diff53prim_5ptrunc_${DATE}.txt

echo "Getting new regions for counts using BID" ${USE_BID}
# get new regions and stats
python ${SRC}/Get_fixed_gene_regions.py --overlap_file ${COUNT_TRUNC_OUT} --new_regions_output_file ${NEW_REG} \
--stats_output_file ${COUNT_STATS} --full_gene_bed ${TRUNC} --use_bids_bedfile ${USE_BID}

wc -l ${NEW_REG}
head -5 ${COUNT_STATS}


###################
# COUNT THE GENES #
###################
echo "Submitted counts with Rsubread for nonoverlapping gene regions by strand......"
#Use this many threads
#Count multi-overlapping read
#Count reads in a strand specific manner
featureCounts \
    -T 32 \
    -O \
    -s 1 \
    -t "exon" \
    -a ${NEW_REG} \
    -F 'GTF' \
    -o ${fixed_counts}/${prefix}_str_fixed_genes_${date}.txt \
    ./*.sorted.bam

wc -l ${fixed_counts}/*

echo "DONE!"
date