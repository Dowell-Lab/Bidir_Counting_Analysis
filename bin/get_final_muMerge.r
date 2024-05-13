library(data.table)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
overlaps_filename = args[1]
tfit_filename = args[2]
dreg_filename = args[3]
output_prefix = args[4]
window = as.numeric(args[5])
sample_ids = args[6]

cat("\nUSING", overlaps_filename, tfit_filename, dreg_filename, output_prefix)
### Reaad in the overlaps between dREG and Tfit
# bedtools intersect -wo -a ${TFIT} -b ${DREG} > overlaps_tfit_dreg.bed
over <- data.frame(fread(overlaps_filename))
dim(over)
#over[1:2,]

# Read in the Tfit & dREG files to determine where the start & overlap columns are
tfit <- fread(tfit_filename)
dreg <- fread(dreg_filename)
# get start & stop columns for dREG
second_start = ncol(tfit)+2
second_start_col = paste0("V", as.character(second_start))
second_stop_col = paste0("V", as.character(second_start+1))
second_sampid_col = paste0("V", as.character(second_start+2))
# get the overlap column
overlap_col = as.character(ncol(tfit) + ncol(dreg) + 1)
## name the regions according to position
over$unique <- paste0(over$V1, ":", over$V2, "-", over$V3)
over$unique_dreg <- paste0(over$V1, ":", over[[second_start_col]], "-", over[[second_stop_col]])
over[1:2,]

###############
## FUNCTIONS ##
###############
# to get the unique set of Sample ids corresponding to a bidirectional if overlapping
get_shrunk <- function(x) {
    x <- strsplit(x, ",")[[1]]
    return(unique(x))
}

###############
## FILTERING ##
###############
cat("\nFILTERING OVERLAPS=============================")

## remove any tfit calls longer than 2.5kb (indicates very low confidence)
over$tfit_len <- over$V3-over$V2
over$dreg_len <- over[[second_stop_col]]-over[[second_start_col]]
over[1:2,]

over <- over[over$tfit_len<2500,]
cat("\nNumber Tfit & dREG overlapping regions saved after removing Tfit calls >2.5kb", length(unique(over$unique)), length(unique(over$unique_dreg)))

## Considering overlaps where > 40% of either section is captured OR MUDIFF < 101
# have unique markers
over$unique_nums <- seq(1, nrow(over))
# get the mus for each
over$tfit_mu <- as.integer((over$V3+over$V2)/2)
over$dreg_mu <- as.integer((over[[second_start_col]]+over[[second_stop_col]])/2)
# get the difference between mus
over$MUDIFF <- abs(over$tfit_mu-over$dreg_mu)
over[1:4,]
cat("\nDifference in mus of overlapping Tfit & dReg calls", 
    quantile(over$MUDIFF, probs = seq (0, 1, 0.25), na.rm = FALSE))
keep <- over[over$MUDIFF < 101,]$unique_nums
cat("\nNumber overlaps kept when only considering overlaps where MUDIFF <101bp", length(keep))

## Remove those with <40% of both regions overlapping
over$tfit_per <- over[[paste0("V", overlap_col)]]/over$tfit_len
over$dreg_per <- over[[paste0("V", overlap_col)]]/over$dreg_len
over_filt <- over[over$tfit_per > 0.4 | over$dreg_per > 0.4,]
cat("\nNumber overlaps kept when removing those with <40% of both regions overlapping", nrow(over_filt))

## Add back regions that <40% both captured BUT MUDIFF < 101
over <- over[over$unique_nums %in% setdiff(keep, over_filt$unique_nums),]
cat("\nNumber of overlaps needing to add back in due to MUDIFF <101bp", nrow(over))
over_filt <- rbind(over_filt, over)

## Get the unique set of sample ids
if (sample_ids) {
    over_filt$SRRs <- paste0(over_filt$V4, ",", over_filt[[second_sampid_col]])
    over_filt$V4 <- as.character(lapply(over_filt$SRRs, get_shrunk))
    }


#######################
## RECOMBINING FILES ##
#######################
cat("\nRECOMBINING TO GET FINAL LIST=============================")

# remove tfit calls with condience intervals larger than 2.5kb
tfit$length <- tfit$V3-tfit$V2
tfit <- tfit[tfit$length < 2500,]

# Remove positions overlapping
tfit$unique <- paste0(tfit$V1, ":", tfit$V2, "-", tfit$V3)
dreg$unique <- paste0(dreg$V1, ":", dreg$V2, "-", dreg$V3)
tfit <- tfit[!tfit$unique %in% over_filt$unique]
dreg <- dreg[!dreg$unique %in% over_filt$unique_dreg]
cat("\nNumber of regions unique to Tfit & dREG", nrow(tfit), nrow(dreg))
cat("\nNumber of regions where overlap", nrow(over_filt))

# recombine with names according to position & algorithm
tfit$name <- paste0(tfit$unique, "-tfit")
dreg$name <- paste0(dreg$unique, "-dreg")
over_filt$name <- paste0(over_filt$unique, "-tfit,dreg")

if (sample_ids) {
    final <- rbind(tfit[,c("V1", "V2", "V3", "name", "V4")], 
               dreg[,c("V1", "V2", "V3", "name", "V4")], 
               over_filt[,c("V1", "V2", "V3", "name", "V4")])
    } else {
    final <- rbind(tfit[,c("V1", "V2", "V3", "name")], 
               dreg[,c("V1", "V2", "V3", "name")], 
               over_filt[,c("V1", "V2", "V3", "name")])
    }

# remove duplicates
final <- final[!duplicated(final[,c("V1", "V2", "V3", "name")])]
# only keep those in primary chromosomes
final <- final[final$V1 %in% c("chr1", "chr2", "chr3", "chr4", "chr5", 
                               "chr6", "chr7", "chr8", "chr9", "chr10", 
                               "chr11", "chr12", "chr13", "chr14", "chr15", 
                               "chr16",  "chr17", "chr18", "chr19", "chr20", 
                               "chr21", "chr22", "chrX", "chrY")]

# recombine with names according to position & algorithm
cat("\nNumber of final regions", nrow(final), length(unique(final$name)))

# # save the file
# write.table(final, paste0(output_prefix, "_MUMERGE_tfit,dreg.bed"), 
#            quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)

# save as according to window
# if no window desired, then just keep the regions as is
if (window==0) {
    write.table(final, paste0(output_prefix, "_MUMERGE_tfit,dreg.bed"), 
            quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
    # if window of 1, then just keep mu
} else if (window==1) {
    final$mu <- as.integer((final$V3+final$V2)/2)
    final$new_start <- final$mu - 1
    if (sample_ids) {
        write.table(final[,c("V1", "new_start", "mu", "name", "V4")], 
           paste0(output_prefix, "_MUMERGE_1bpwin_tfit,dreg.bed"), 
           quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
        } else {
        write.table(final[,c("V1", "new_start", "mu", "name")], 
           paste0(output_prefix, "_MUMERGE_1bpwin_tfit,dreg.bed"), 
           quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
        }
    
} else {
    window_half = as.integer(window/2)
    final$mu <- as.integer((final$V3+final$V2)/2)
    final$Vb <- final$mu + window_half
    final$Va <- final$mu - window_half
    # if any are below 1 position, fix to be at 1
    final[final$Va < 1,]$Va <- 1
    if (sample_ids) {
        write.table(final[,c("V1", "Va", "Vb", "name", "V4")], 
           paste0(output_prefix, "_MUMERGE_", window, "bpwin_tfit,dreg.bed"), 
           quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
        } else {
        write.table(final[,c("V1", "Va", "Vb", "name")], 
           paste0(output_prefix, "_MUMERGE_", window, "bpwin_tfit,dreg.bed"), 
           quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
        }
    
    }





