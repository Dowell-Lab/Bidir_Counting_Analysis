library(data.table)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
mumerge_filename = args[1]
count_win = as.numeric(args[2])
tss_win = as.numeric(args[3])
output_prefix = args[4]

cat("\nUSING File", mumerge_filename,"\nCount & TSS windows", count_win, tss_win, "\nOutput Prefix", output_prefix)
### Read in the consensus region
# bedtools intersect -wo -a ${TFIT} -b ${DREG} > overlaps_tfit_dreg.bed
cons <- data.frame(fread(mumerge_filename))
dim(cons)
#over[1:2,]

# Get the names for each bidir
if (nrow(cons) > 3) {
    cons$unique_name = cons$V4
    } else {
    cons$unique_name = paste0(cons$V1, ":", cons$V2, "-", cons$V3) 
    }


###############
## FILTERING ##
###############
cat("\nFILTERING Poor Confidence Regions=============================")

## remove any tfit calls longer than 3.5kb (indicates very low confidence)
cons$length <- cons$V3-cons$V2
cons <- cons[cons$length<3500,]
cat("\nNumber regions saved after removing calls >3.5kb -- not trustworthy", length(unique(cons$unique)), length(unique(cons$unique_dreg)))

###############
## Save according to windows ##
###############
# get mu
cons$mu <- as.integer((cons$V2+cons$V3)/2)
cons$count_start <- cons$mu - count_win
cons$count_stop <- cons$mu + count_win
cons$tss_start <- cons$mu - tss_win
cons$tss_stop <- cons$mu + tss_win

write.table(cons[,c("V1", "count_start", "count_stop", "unique_name")], 
                paste0(output_prefix, "_MUMERGE_", count_win, "win_count.sorted.bed"), 
            quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
write.table(cons[,c("V1", "tss_start", "tss_stop", "unique_name")], 
                paste0(output_prefix, "_MUMERGE_", count_win, "win_TSS.sorted.bed"), 
            quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)







