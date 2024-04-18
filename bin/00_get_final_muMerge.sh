#!/bin/bash

source activate sclinker
bedtools --version
##########################
# EDIT THE FOLLOWING
##########################
# files
INPUT_DIR=/Users/hopekirby/Desktop/Resp_Genetics/Asthma_Genetics/Sasse2023/beds
OUTPUT_DIR=/Users/hopekirby/Desktop/TF_Enh_Linking/eRNA-SNP/out
TFIT=${INPUT_DIR}/Sasse2019nascent_hg38_tfit_MUMERGE.bed
DREG=${INPUT_DIR}/Sasse2019nascent_hg38_dreg_MUMERGE.bed
# prefix
OUT_PREFIX=Sasse2019nascent_04-4-24
# TODO: allow multiple windows to be done at the same time
# fixed window you want to have (0= just use confidence intervals, 1=just mu, 300= 150± mu)
WINDOW=50


##########################
# STARTING THE CODE
##########################
# make an overlaps directory if it doesn't already exist
mkdir -p ${OUTPUT_DIR}/overlaps
OVER_OUT=${OUTPUT_DIR}/overlaps/overlaps_${OUT_PREFIX}.bed
echo "====== Getting overlaps between Tfit & dREG: ", ${TFIT}, ${DREG}
# Get the overlaps between the Tfit muMerge
bedtools intersect -wo -a ${TFIT} -b ${DREG} > ${OVER_OUT}
head ${OVER_OUT}
echo "Saved at" ${OVER_OUT}

echo "====== Getting the final muMerge file"
## arguments: overlaps_filename, tfit_filename, dREG_filename, output_prefix, dREG start col, overlap col 
# move output to the beds folder
OUT_PREFIX=${OUTPUT_DIR}/beds/${OUT_PREFIX}
Rscript /Users/hopekirby/Desktop/TF_Enh_Linking/eRNA-SNP/00_Prelim/get_final_muMerge.r ${OVER_OUT} ${TFIT} ${DREG} ${OUT_PREFIX} ${WINDOW}
echo "\nFinal files saved with prefix" ${OUT_PREFIX}
ls -lh ${OUT_PREFIX}*

## perform bedtools sort
bedtools sort -i ${OUT_PREFIX}_MUMERGE_${WINDOW}bpwin_tfit,dreg.bed > ${OUT_PREFIX}_MUMERGE_${WINDOW}bpwin_tfit,dreg.sorted.bed
wc -l ${OUT_PREFIX}_MUMERGE_${WINDOW}bpwin_tfit,dreg.bed
wc -l ${OUT_PREFIX}_MUMERGE_${WINDOW}bpwin_tfit,dreg.sorted.bed
rm ${OUT_PREFIX}_MUMERGE_${WINDOW}bpwin_tfit,dreg.bed


echo "DONE!"
date