
##########################
# EDIT THE FOLLOWING
##########################
source activate (conda environment with R, R data.table, and R stringr)
# working directory (overlaps in wd/overlaps, regions saved in wd/regions)
WD=
# date and prefix ussed previously for the overlaps
PREFIX=
DATE=
#original bidirectionals for counting (reccomend 600bp window for now)
OG_BID=
# minimum rowSUM counts needed for a gene to be considered "transcribed" and significantly altering bid counts (I used 40 for 4 samples)
COUNT_LIMIT=40

##########################
# STARTING THE CODE
##########################
mkdir -p ${WD}/regions
Rscript get_TSS_and_filt_bids.r ${WD} ${PREFIX} ${DATE} ${OG_BID} ${COUNT_LIMIT}
ls ${WD}/regions