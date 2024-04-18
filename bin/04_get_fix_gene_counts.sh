#!/bin/bash

#SBATCH --job-name=featureCounts_bids_SE
#SBATCH --output=/scratch/Users/hoto7260/Resp_Env/Sasse2019/e_and_o/featurecounts_bids_%j.out
#SBATCH --error=/scratch/Users/hoto7260/Resp_Env/Sasse2019/e_and_o/featurecounts_bids_%j.err
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=5G
#SBATCH --partition short
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hoto7260@colorado.edu

##load modules
module load samtools/1.8
module load subread/1.6.2

source activate /Users/hoto7260/miniconda3/envs/python39 (Conda environment with python 3.9, pandas (I used 2.0.3), numpy (I used 1.25.2)

# PRINT THE VERSIONS USED
python --version
conda list

##########################
# EDIT THE FOLLOWING
##########################
## Naming
prefix="Sasse2019_nascent"
date=04_04_24
## IN/OUT Directories
SRC= # where the github repo is
# wd is where the bams and counts will be saved
WD=/scratch/Users/hoto7260/Resp_Env/Sasse2019/out

TRUNC="/scratch/Shares/dowell/genomes/hg38/ncbi/hg38_refseq_diff53prime_5ptrunc_counting.bed"
TRUNC_OVERLAPS="/scratch/Users/hoto7260/Resp_Env/Sasse2019/out/overlaps/overlaps_hg38_trunc_Sasse2019nascent_MUMERGE_600bp_04-4-24.bed"
USE_BID=${wd}/regions/${prefix}_nontssbid_uniqueid_above25_${date}.bed


###################################
# GET THE FIXED REGIONS FOR GENES #
###################################
bams=${WD}/bams
counts=${WD}/counts
fixed_counts=${WD}/fixed_counts 
NEW_REG=${wd}/regions/${prefix}_new_regions_len_posneg_gene_filt_divconv_diff53prim_5ptrunc_${date}.gtf
COUNT_STATS=${wd}/regions/${prefix}_stats_len_posneg_gene_filt_divconv_diff53prim_5ptrunc_${date}.txt

echo "Getting new regions for counts using BID" ${USE_BID}
# get new regions and stats
python ${SRC}/bin/Get_fixed_gene_regions.py --overlap_file ${TRUNC_OVERLAPS} --new_regions_output_file ${NEW_REG} \
--stats_output_file ${COUNT_STATS} --full_gene_bed ${TRUNC} --use_bids_bedfile ${USE_BID}

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
    ${bams}/*.sorted.bam
