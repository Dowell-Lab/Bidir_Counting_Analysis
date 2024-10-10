#!/bin/bash

#SBATCH --job-name=Bidir_Counting
#SBATCH --output=/scratch/Users/path/e_and_o/comb_get_counts_%j.out
#SBATCH --error=/scratch/Users/path/e_and_o/comb_get_counts_%j.err
#SBATCH --time=03:00:00
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
WD=/scratch/Users/
### Naming
PREFIX=
DATE=

## IN Directories Need mmfiltbams. If mmfiltbams is empty, must have crams or bams
##   NOTE: MMFILTbams are usually produced by Bidirectional-Flow pipeline
##   NOTE: The mmfiltbams directory can have bams that are not multimapped filtered, as long as the bams of interest have the suffix mmfilt.sorted.bam
mmfiltbams=""
crams=""
bams=""

### Files
# File with consensus regions (e.g. output from Mumerge)
CONS_FILE=

### Parameters
# Type of counting (MU_COUNTS, SIMPLE, or BOTH)
COUNT_TYPE="MU_COUNTS"
# window for counting (mu - window & mu + window)
COUNT_WIN=1000
# window for TSS overlaps (mu - window & mu + window)
TSS_WIN=25
# minimum coverage needed for a gene to be considered "transcribed" and significantly altering bid counts 
COUNT_LIMIT_GENES=0.7
# the number required for a bidirectional to be considered impeding on gene transcription
# usually do 25 for 4 samples
COUNT_LIMIT_BIDS=30
# do you want to get TFEA combined file outputs? YES or NO
TFEA="YES"
# do you want to have unique names be created for each bidirectional?
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
echo "min coverage for gene to be a convolution problem" ${COUNT_LIMIT_GENES}
echo "min rowSum counts for a bid to be a convolution problem" ${COUNT_LIMIT_BIDS}
echo ""

echo "##########################################################################"
echo "##################         Getting Bams ready         ####################"
echo "##########################################################################"
# GETTING THE BAMS IF NEEDED
# if there is no mmfiltbams
if [[ -z "$mmfiltbams" ]]; then
    if [[ -n "$crams" ]]; then
        echo "Getting multimapped filtered bams from crams"
        mmfiltbams=${WD}/mapped/mmfiltbams
        mkdir -p ${mmfiltbams}
        for cram in ${crams}/*.sorted.cram; do
            prefix=$(basename "$cram" | sed 's/\.sorted\.cram$//')
            samtools view -@ 32 -b -1 -T ${genome} ${cram} | samtools view -@ 32 -h -q 1 | \
               grep -P '(NH:i:1|^@)' | \
               samtools view -h -b > ${mmfiltbams}/${prefix}.mmfilt.sorted.bam
            samtools index ${mmfiltbams}/${prefix}.mmfilt.sorted.bam ${mmfiltbams}/${prefix}.mmfilt.sorted.bam.bai
        done
    elif [[ -n "$bams" ]]; then
        echo "Getting multi-mapped filtered bams from bams"
        mmfiltbams=${WD}/mapped/mmfiltbams
        mkdir -p ${mmfiltbams}
        for bam in ${bams}/*.sorted.bam; do
            prefix=$(basename "$bam" | sed 's/\.sorted\.bam$//')
            samtools view -@ 32 -h -q 1 ${bam} | \
               grep -P '(NH:i:1|^@)' | \
               samtools view -h -b > ${mmfiltbams}/${prefix}.mmfilt.sorted.bam
           samtools index ${mmfiltbams}/${prefix}.mmfilt.sorted.bam ${mmfiltbams}/${prefix}.mmfilt.sorted.bam.bai
        done
    else; then
        echo "Need to have at least one of the variables filled: mmfiltbams, crams, bams"
        exit 152
    fi
else; then
    echo "Using the provided multimapped filtered bams at ${mmfiltbams}"
fi


echo "#####################################################################################"
echo "########################            1a. GET THE OVERLAPS         ####################"
echo "#####################################################################################"

# The FULL transcripts with putatives (WARNING: if names are too long, bedtools breaks without warning). This has been addressed with the _fixnames file.
TRANSCRIPTS=${SRC}/assets/hg38_refseq_diff53prime_with_putatives_fixnames.sorted.bed
# The gene file you will use for counting (make sure truncated)
TRUNC=${SRC}/assets/hg38_refseq_diff53prime_5ptrunc_counting.sorted.bed
# the gene file where there are 
GENE_DWN=${SRC}/assets/hg38_refseq_diff53prime_10kbdwnstm_with_putatives_fixnames.sorted.bed
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
DWN_OUT=${OUT_DIR}/overlaps_hg38_withput_dwnstm_${PREFIX}_MUMERGE_${DATE}.bed
COUNT_TRUNC_OUT=${OUT_DIR}/overlaps_hg38_trunc_${PREFIX}_MUMERGE_${DATE}.bed
CLOSE_OUT=${OUT_DIR}/closest_hg38_withput_${PREFIX}_MUMERGE_${DATE}.bed
CLOSE_BID=${OUT_DIR}/closest_Bids_${PREFIX}_MUMERGE_${DATE}.bed

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
## Bids that overlap with transcripts -- must have gene first 
bedtools intersect -wo -a ${TRANSCRIPTS} -b ${COUNT_WIN_FILE} > ${COUNT_OUT}
wc -l ${COUNT_OUT}
## Bids that overlap with 10kb downstream of transcripts -- must have gene first 
bedtools intersect -wo -a ${GENE_DWN} -b ${COUNT_WIN_FILE} > ${DWN_OUT}
wc -l ${DWN_OUT}
## Bids that overlap with truncated transcripts -- must have gene first
bedtools intersect -wo -a ${TRUNC} -b ${COUNT_WIN_FILE} > ${COUNT_TRUNC_OUT}
wc -l ${COUNT_TRUNC_OUT}
## Bids possibly overlapping with one another (-N means can't have same name)
bedtools closest -D ref -k 5 -N -a ${TSS_WIN_FILE} -b ${TSS_WIN_FILE} > ${CLOSE_BID}

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


echo "Running bedtools coverage to see the coverage over a gene with multimap filtered bams"
# use bedtools coverage to ensure that most of the gene is being transcribed
cd ${mmfiltbams}
bedtools coverage -s -sorted -a ${TRUNC} -b *mmfilt.sorted.bam > ${counts}/${PREFIX}_str_put_genes.txt
head -3 ${counts}/${PREFIX}_str_put_genes.txt

    
echo ""
echo "###################################################################################"
echo "######     2. GET TSS Bids & GTF files for Bid Counting (no overlaps)       #######"
echo "###################################################################################"
# YES refers to TFEA output being mades
Rscript ${SRC}/bin/get_TSS_and_filt_bids.r ${WD} ${PREFIX} ${DATE} ${COUNT_LIMIT_GENES} ${COUNT_WIN} ${COUNT_TYPE} ${TFEA}
wc -l ${WD}/regions/*


echo ""
echo "###################################################################################"
echo "########################         3. GET Fixed & Bid Counts     ####################"
echo "###################################################################################"
cd ${mmfiltbams}

if [[ "$COUNT_TYPE" == "SIMPLE" || "$COUNT_TYPE" == "BOTH" ]]; then
    bidirs=${WD}/regions/${PREFIX}_MUMERGE_tfit,dreg_${DATE}_filt.saf
    echo "Submitting counts for Bidirectionals with Simple SAF"
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
       ./*mmfilt.sorted.bam
    
    # Count reads NOT in a strand specific manner
    featureCounts \
       -T 32 \
       -O \
       -s 1 \
       -a ${bidirs} \
       -F 'SAF' \
       -o ${counts}/${PREFIX}_pos_bidirs.txt \
       ./*mmfilt.sorted.bam
    
    # Count multi-overlapping read
    featureCounts \
       -T 32 \
       -O \
       -s 2 \
       -a ${bidirs} \
       -F 'SAF' \
       -o ${counts}/${PREFIX}_neg_bidirs.txt \
       ./*mmfilt.sorted.bam

    # Fix the counts
    Rscript ${SRC}/bin/fix_bid_counts.r ${WD} ${PREFIX} ${DATE} ${COUNT_LIMIT_BIDS} ${COUNT_LIMIT_GENES} "SIMPLE"
fi

if [[ "$COUNT_TYPE" == "MU_COUNTS" || "$COUNT_TYPE" == "BOTH" ]]; then
    gtf_prefix=${WD}/regions/${PREFIX}_mucounts_${COUNT_WIN}_${DATE}
    # count using the GTFs
    featureCounts \
       -T 32 \
       -O \
       -s 1 \
       -a ${gtf_prefix}.gtf \
       -t "exon" \
       -F 'GTF' \
       -o ${counts}/${PREFIX}_mucounts_str_bidirs.txt \
       ./*mmfilt.sorted.bam

    featureCounts \
       -T 32 \
       -O \
       -s 1 \
       -a ${gtf_prefix}_pos.gtf \
       -t "exon" \
       -F 'GTF' \
       -o ${counts}/${PREFIX}_mucounts_pos_bidirs.txt \
       ./*mmfilt.sorted.bam

       featureCounts \
       -T 32 \
       -O \
       -s 1 \
       -a ${gtf_prefix}_neg.gtf \
       -t "exon" \
       -F 'GTF' \
       -o ${counts}/${PREFIX}_mucounts_neg_bidirs.txt \
       ./*mmfilt.sorted.bam

       # Fix the counts
       Rscript ${SRC}/bin/fix_bid_counts.r ${WD} ${PREFIX} ${DATE} ${COUNT_LIMIT_BIDS} ${COUNT_LIMIT_GENES} "MU_COUNTS"
fi


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
python ${SRC}/bin/Get_fixed_gene_regions.py --overlap_file ${COUNT_TRUNC_OUT} --new_regions_output_file ${NEW_REG} \
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
cd ${mmfiltbams}
featureCounts \
    -T 32 \
    -O \
    -s 1 \
    -t "exon" \
    -a ${NEW_REG} \
    -F 'GTF' \
    -o ${fixed_counts}/${prefix}_str_fixed_genes_${date}.txt \
   ./*mmfilt.sorted.bam

wc -l ${fixed_counts}/*

###################
# REMOVE INTERMEDIATE FILES #
###################
rm ${TSS_WIN_FILE}
rm ${COUNT_WIN_FILE}

echo "DONE!"
date
