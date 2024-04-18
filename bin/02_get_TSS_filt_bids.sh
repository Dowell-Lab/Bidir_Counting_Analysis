
##########################
# EDIT THE FOLLOWING
##########################
source activate (conda environment with R, R data.table, and R stringr)
# directory wherre overlaps are
OVER_DIR=
# date and prefix ussed previously for the overlaps
PREFIX=
DATE=
#original bidirectionals for counting (reccomend 600bp window for now)
OG_BID=
# directory where save TSS bids & SAF for bid counting
OUT_DIR=
# minimum rowSUM counts needed for a gene to be considered "transcribed" and significantly altering bid counts (I used 40 for 4 samples)
COUNT_LIMIT=40


##########################
# STARTING THE CODE
##########################
Rscript get_TSS_and_filt_bids.r ${OVER_DIR} ${PREFIX} ${DATE} ${OG_BID} ${OUT_DIR} ${COUNT_LIMIT}
ls ${OUT_DIR}