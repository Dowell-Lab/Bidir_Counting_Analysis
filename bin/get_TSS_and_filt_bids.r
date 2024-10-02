library(data.table)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
cat("\nArguments: ", args, "\n")
wd = args[1] # working directory
prefix = args[2] #prefix for files
date = args[3] #
count_limit = as.integer(args[4])
count_win = as.integer(args[5])
type = args[6] # type of regions file for counting (SIMPLE or MU_COUNTS or BOTH)
TFEA = args[7]

og_bid_file=paste0(wd, "/regions/", prefix, "_", date, "_MUMERGE_", count_win, "win_count.sorted.bed")
out_dir = paste0(wd, "/regions/")
overlaps_dir = paste0(wd, "/overlaps/")
tss_out_file = paste0(overlaps_dir, "overlaps_hg38_TSS1kb_withput_", prefix, "_MUMERGE_", date, ".bed")
count_out_file = paste0(overlaps_dir, "overlaps_hg38_withput_", prefix, "_MUMERGE_", date, ".bed")
closest_file = paste0(overlaps_dir, "closest_hg38_withput_", prefix, "_MUMERGE_", date, ".bed")
count_file = paste0(wd, "/counts/", prefix, "_str_put_genes.txt")
close_bid_file = paste0(overlaps_dir, "closest_Bids_", prefix, "_MUMERGE_" date, ".bed")


###########################
#####    FUNCTIONS    #####
###########################

get_info_for_calls <- function(close_df, fixed_length=1000) {
    # This function gets the information needed in the overlaps file to get the
    #     appropriate leftmost and rightmost calls
    # PARAMETERS:
    #     close_df: a dataframe with the columns chrom, start, stop, name of region 1, the same for region 2, and finally the distance between the two regions. It is assumed that the middle of the start and stop of each region (mu) is the starting point of each transcript.
    #     fixed_length: Fixed length distance from mu to consider the "original" length of the transcripts. This fixed_length + 100 is also the max distance the mus of two features can be from one another to consider potential overlap (since the count region wouldn't get farther anyways)
    diff_call = fixed_length+100
    close_df$mu1 <- as.integer((close_df$V2+close_df$V3)/2)
    close_df$mu2 <- as.integer((close_df$V6+close_df$V7)/2)
    close_df$MUDIFF <- close_df$mu1-close_df$mu2
    close_df <- close_df[abs(close_df$MUDIFF) < diff_call,]
    # get the desired distances of interest 
    close_df$V2 = close_df$mu1 - fixed_length
    close_df$V3 = close_df$mu1 + fixed_length
    close_df$V6 = close_df$mu2 - fixed_length
    close_df$V7 = close_df$mu2 + fixed_length
    
    # get leftmost and rightmost
    close_df$leftmost_Bidir <- "NA"; close_df$leftmost_left <- 0; close_df$leftmost_mu <- 0; close_df$leftmost_right <- 0
    close_df$rightmost_Bidir <- "NA"; close_df$rightmost_right <- 0; close_df$rightmost_mu <- 0; close_df$rightmost_left <- 0
    # where first is leftmost one
    close_df[close_df$MUDIFF < 0,]$leftmost_Bidir <- close_df[close_df$MUDIFF < 0,]$V4
    close_df[close_df$MUDIFF < 0,]$leftmost_left <- close_df[close_df$MUDIFF < 0,]$V2
    close_df[close_df$MUDIFF < 0,]$leftmost_mu <- close_df[close_df$MUDIFF < 0,]$mu1
    close_df[close_df$MUDIFF < 0,]$rightmost_Bidir <- close_df[close_df$MUDIFF < 0,]$V8
    close_df[close_df$MUDIFF < 0,]$rightmost_right <- close_df[close_df$MUDIFF < 0,]$V7
    close_df[close_df$MUDIFF < 0,]$rightmost_mu <- close_df[close_df$MUDIFF < 0,]$mu2
    close_df[close_df$MUDIFF < 0,]$leftmost_right <- close_df[close_df$MUDIFF < 0,]$V3
    close_df[close_df$MUDIFF < 0,]$rightmost_left <- close_df[close_df$MUDIFF < 0,]$V6
    # where second is leftmost one
    close_df[close_df$MUDIFF > 0,]$rightmost_Bidir <- close_df[close_df$MUDIFF > 0,]$V4
    close_df[close_df$MUDIFF > 0,]$rightmost_right <- close_df[close_df$MUDIFF > 0,]$V3
    close_df[close_df$MUDIFF > 0,]$rightmost_mu <- close_df[close_df$MUDIFF > 0,]$mu1
    close_df[close_df$MUDIFF > 0,]$leftmost_Bidir <- close_df[close_df$MUDIFF > 0,]$V8
    close_df[close_df$MUDIFF > 0,]$leftmost_left <- close_df[close_df$MUDIFF > 0,]$V6
    close_df[close_df$MUDIFF > 0,]$leftmost_mu <- close_df[close_df$MUDIFF > 0,]$mu2
    close_df[close_df$MUDIFF > 0,]$leftmost_right <- close_df[close_df$MUDIFF > 0,]$V7
    close_df[close_df$MUDIFF > 0,]$rightmost_left <- close_df[close_df$MUDIFF > 0,]$V2
    # check by new MUDIFF should all be positive
    close_df$leftmost_mu <- as.integer(close_df$leftmost_mu)
    close_df$rightmost_mu <- as.integer(close_df$rightmost_mu)
    close_df$newMUDIFF <- close_df$rightmost_mu - close_df$leftmost_mu
    if (nrow(close_df[close_df$newMUDIFF < 0,]) > 0) {
        stop("There is a problem in the overlap files of the Bids") }
    
    return(close_df)
    }

get_new_pos <- function(close_df) {
    # This function gets the new positions for regions BUT assumes the regions are only overlapping
    #      one other region if any.
    # PARAMETERS
    #     close_df: a dataframe with the columns chrom, start, stop, name of region 1, the same for region 2, and finally the distance between the two regions. It is assumed that the middle of the start and stop of each region (mu) is the starting point of each transcript.
    # make easier naming of info
    close_df$left_neg_stop = close_df$leftmost_mu 
    close_df$left_neg_start = close_df$leftmost_left
    close_df$right_neg_start = pmax(close_df$leftmost_mu, close_df$rightmost_left) # new right neg 'start' should be rightmost of left mu or the original end
    close_df$right_neg_stop = close_df$rightmost_mu
    close_df$left_pos_start = close_df$leftmost_mu
    close_df$left_pos_stop = pmin(close_df$rightmost_mu, close_df$leftmost_right) # new right pos should end at min of right mu or original call
    close_df$right_pos_start = close_df$rightmost_mu
    close_df$right_pos_stop = close_df$rightmost_right
    print(dim(close_df[(close_df$left_neg_stop - close_df$left_neg_start <= 0) | 
                       (close_df$left_pos_stop - close_df$left_pos_start <= 0) | 
                       (close_df$right_neg_stop - close_df$right_neg_start <= 0) | 
                       (close_df$right_pos_stop - close_df$right_pos_start <= 0),]))
    return(close_df)
    }

get_inbetween_lines <- function(close_df, in_between) {
    # This function gets the GTF lines for cases where a region has overlapping regions
    #     on both sides of it.
    # PARAMETERS
    #     close_df: a dataframe with the columns chrom, start, stop, name of region 1, the same for region 2, and finally the distance between the two regions. It is assumed that the middle of the start and stop of each region (mu) is the starting point of each transcript.
    #     in_between: the list of regions that have overlapping regions on both sides 
    # for each in between
    between_lines = c()
    for (between_bid in in_between) {
        # get the features closest where BID is on the left
        left = close_df[close_df$leftmost_Bid == between_bid,]
        # get the features closest where BID is on the right
        right = close_df[close_df$rightmost_Bid == between_bid,]
        # get the gene coord (just left and right mu)
        # get the pos coord (mu & leftest right mu on right)
        # get the neg coord (rightest left mu & mu)
        between_lines = c(between_lines, paste0(left$V1[1], "\t.\tgene\t", left$leftmost_left[1], 
                                                "\t", left$leftmost_right[1], '\t.\t.\t.\tgene_id "', between_bid, '";'), 
                                              paste0(left$V1[1], "\t.\texon\t", left$leftmost_mu[1], "\t",
                                              min(c(left$rightmost_mu, left$leftmost_right[1])),  '\t.\t+\t.\tgene_id "', between_bid, '"; transcript_id "', between_bid, '";'), 
                                              paste0(left$V1[1], "\t.\texon\t", max(c(right$leftmost_mu[1], left$leftmost_left[1])), "\t", right$rightmost_mu[1],  '\t.\t-\t.\tgene_id "', between_bid, 
                                                  '"; transcript_id "', between_bid, '";'))
    }
    return(between_lines)
    }



get_GTF_lines <- function(close_df) {
    # This function gets the annotation lines to use for a GTF that accounts for overlapping regions, including regions that have overlap with other regions on both strands.
    # PARAMETERS
    #       close_df: close_df that has already gone through the functions get_info_for_calls and get_new_pos.
    in_between <- intersect(close_df$rightmost_Bidir, close_df$leftmost_Bidir)
    cat("\nIn between", length(in_between), "\n")
    inbetween_lines = get_inbetween_lines(close_df, in_between)
    close_df$left_gene_row = paste0(close_df$V1, "\t.\tgene\t", close_df$leftmost_left, "\t", close_df$leftmost_right, 
                                   '\t.\t.\t.\tgene_id "', close_df$leftmost_Bidir, '";')
    close_df$left_neg_row = paste0(close_df$V1, "\t.\texon\t", close_df$left_neg_start, "\t", close_df$left_neg_stop, 
                                   '\t.\t-\t.\tgene_id "', close_df$leftmost_Bidir, 
                                  '"; transcript_id "', close_df$leftmost_Bidir, '";' )
    close_df$left_pos_row = paste0(close_df$V1, "\t.\texon\t", close_df$left_pos_start, "\t", close_df$left_pos_stop, 
                                   '\t.\t+\t.\tgene_id "', close_df$leftmost_Bidir, 
                                  '"; transcript_id "', close_df$leftmost_Bidir, '";' )
    close_df$right_gene_row = paste0(close_df$V1, "\t.\tgene\t", close_df$rightmost_left, "\t", close_df$rightmost_right, 
                                   '\t.\t.\t.\tgene_id "', close_df$rightmost_Bidir, '";')
    close_df$right_neg_row = paste0(close_df$V1, "\t.\texon\t", close_df$right_neg_start, "\t", close_df$right_neg_stop, 
                                   '\t.\t-\t.\tgene_id "', close_df$rightmost_Bidir, 
                                  '"; transcript_id "', close_df$rightmost_Bidir, '";' )
    close_df$right_pos_row = paste0(close_df$V1, "\t.\texon\t", close_df$right_pos_start, "\t", close_df$right_pos_stop, 
                                   '\t.\t+\t.\tgene_id "', close_df$rightmost_Bidir, 
                                  '"; transcript_id "', close_df$rightmost_Bidir, '";' )
    df <- close_df[,c("left_gene_row", "left_neg_row", "left_pos_row", 
                 "right_gene_row", "right_neg_row", "right_pos_row")]
    cat("\nNumber of rows so far", dim(df))
    lines_vector <- as.vector(t(df))
    #between_removing_regexp <- paste(in_between, collapse = "|")
    #return(list("bet"=between_removing_regexp, "lines"=lines_vector))
    cat("\nNumber of lines so far", length(lines_vector))
    # if some in between
    if (length(in_between) != 0) {
        # the in_between into chunks of 100
        chunks <- split(my_vect, ceiling(seq_along(my_vect) / 100))
        for (chunk in chunks) {
             between_removing_regexp <- paste(chunk, collapse = "|")
            lines_vector <- lines_vector[!grepl(between_removing_regexp, lines_vector)]
            }
        # add the correct lines for the in between
        lines_vector <- c(lines_vector, inbetween_lines)
        }
    
    return(lines_vector)
    }


###################################s
## 1. Get the TSS Bidirectionals ##
###################################
# read in the TSS region overlaps
### NOTE: For now, This must be in the format so the gene name is at V8
over <- fread(tss_out_file)
# Get the distance between mu and TSS
over$mu <- as.integer((over$V3+over$V2)/2)
over$tss <- as.integer((over$V6+over$V7)/2)
over$MUDIFF <- abs(over$mu-over$tss)

## Get the Gene (not transcript) name
over$Gene <- str_split_fixed(over$V8, ":", 2)[,1]
trans_with_tss <- unique(over$V8)
genes_with_tss <- unique(over$Gene)

## Get the TSS bidirectionals (note if hand annotation might be needed)
over <- data.frame(over)
over$unique <- seq(1, nrow(over))
hand_annotate = c(); TSS_used = c(); keep <- c(); remove <- c()
for (isoform in trans_with_tss) {
    # get the overlaps
    filt = over[over$V8 == isoform,]
    # get the TSS with the minimum MUDIFF
    tss = filt[filt$MUDIFF == min(filt$MUDIFF),]$V4
    #cat("\n", isoform, tss)
    # if more than one bidirectional has the same min distance --> hand annotate (for now just keep first)
    if (length(tss) != 1) {
        hand_annotate <- c(hand_annotate, isoform)
        cat("\nThe isoform", isoform, "has", length(tss), "TSS that are closest. Hand annotate")
        tss = tss[1]
    } 
    if (tss %in% TSS_used) {
        # if the TSS has already been assigned to a gene
        # check the difference in MUDIFF between the two genes
        filt2 = over[over$V4 == tss,]
        # if the MUDIFF of the two genes are significantly different, only keep one
        if (max(abs(filt2$MUDIFF))-min(abs(filt2$MUDIFF)) > 50) {
            keep <- c(keep, filt2[abs(filt2$MUDIFF) == min(abs(filt2$MUDIFF)),]$unique)
            remove <- c(remove, filt2[abs(filt2$MUDIFF) == max(abs(filt2$MUDIFF)),]$unique)
        } else {
            keep <- c(keep, filt$unique)
        }}else {
            keep <- c(keep, filt[filt$V4 == tss,]$unique)
            }
    TSS_used <- union(TSS_used, tss)
        }
# only address the TSS bidirectionals clarified to keep and remove
over_filt <- over[over$unique %in% keep,]
over_filt <- over_filt[!over_filt$unique %in% remove,]
# Record the numbers changed
trans_with_tss_filt <- unique(over_filt$V8)
cat("\nOriginal # Transcripts captured:", length(trans_with_tss), "and new #:", length(trans_with_tss_filt))
genes_with_tss_filt <- unique(over_filt$Gene)
cat("\nOriginal # Genes captured:", length(genes_with_tss), "and new #:", length(genes_with_tss_filt))
cat("\nOriginal MUDIFF quantiles vs Filtered:")
quantile(over$MUDIFF, probs = seq (0, 1, 0.25), na.rm = FALSE)
quantile(over_filt$MUDIFF, probs = seq (0, 1, 0.25), na.rm = FALSE)

# Save the bidirectionals
colnames(over_filt) <- c("Chr", "Bid_Start", "Bid_Stop", "BidID", "Gene_Chr", "Gene_Start", "Gene_Stop", 
                        "TranscriptID", ".", "Strand", "Length", "mu", "TSS", "MUDIFF", "GeneID", "unique")
over_filt <- over_filt[,c("Chr", "Bid_Start", "Bid_Stop", "BidID", "Gene_Start", "Gene_Stop",
             "TranscriptID", "Strand", "mu", "TSS", "GeneID")]
over_filt[1:2,]

write.table(over_filt, paste0(out_dir, "tss_bid_", prefix, "_", date, ".txt"), 
           row.names=FALSE, sep="\t", quote=FALSE)

##############################################
## 2. Filter the bids according to overlaps ##
##############################################
# read in the overlaps for counting
count_over <- fread(count_out_file)
# read in the prel gene counts
counts <- fread(count_file)
# ensure that the Genes and overlaps have the same Geneid format
counts$Geneid <- str_split_fixed(counts$Geneid, ",", 2)[,1]
length(setdiff(count_over$V4, counts$Geneid))
# get the low counts
count_matrix = data.frame(counts[,7:ncol(counts)])
low_counts = rownames(count_matrix[rowSums(count_matrix) < count_limit,])
# only consider overlaps with genes with OUT low counts
count_over <- count_over[!count_over$V4 %in% low_counts]

# Split the overlaps into positive and negative strands
count_over_neg = count_over[count_over$V6 == "-", ]
count_over_pos = count_over[count_over$V6 == "+", ]

pos_bids <- unique(count_over_pos$V10)
neg_bids <- unique(count_over_neg$V10)
pos_neg_filt <- intersect(pos_bids, neg_bids)
# remove those with TSS
pos_neg_filt <- setdiff(pos_neg_filt, over_filt$BidID)
cat("\nNumber of NonTSS bids with overlap of transcribed genes on both strands", length(pos_neg_filt))

#############################################################################
## 3. Filter the bids according to 5kb downstream of TWO transcribed genes ##
#############################################################################
closest <- fread(closest_file)
# only keep those with 10kb distance from either end
closest <- closest[abs(closest$V11) < 10000,]
# only keep those with genes that are being transcribed significantly
closest <- closest[!closest$V4 %in% low_counts]
# only keep those that are within 10kb to the 3prime end 
# for positive genes, means V11 <10000 & positive 
closest <- closest[(closest$V6 == "+" & closest$V11 > 0 & closest$V11 < 10001) | closest$V6 == "-",]
nrow(closest)
# for negative genes, means V11 <10000 & negative
closest <- closest[(closest$V6 == "-" & closest$V11 < 0 & closest$V11 > -10001) | closest$V6 == "+",]
# now see which bids are this for two colliding genes and NOT TSS
pos_colliding <- closest[closest$V6 == "+",]$V10
neg_colliding <- closest[closest$V6 == "-",]$V10
pos_neg_colliding <- setdiff(intersect(pos_colliding, neg_colliding), 
                             over_filt$BidID)
cat("\nNumber of NonTSS bids within 5kb of 3' of transcribed genes on both strands", length(pos_neg_colliding))

# read in the original bids
bids <- fread(og_bid_file)
pos_neg_remove <- union(pos_neg_colliding, pos_neg_filt)
cat("\nTotal Number of NonTSS bids removing due to convolution", length(pos_neg_remove))
bids <- bids[!bids$V4 %in% pos_neg_remove]
cat("\nTotal Number of Bids remaining", nrow(bids))


######################################################
## 4. Get the counting files for the remaining bids ##
######################################################
### SIMPLE SAF ACCORDING TO FIXED WINDOW LENGTH
if (type == "SIMPLE" | type == "BOTH") {
    # save the remaining bids as a saf file
    colnames(bids) <- c("Chr", "Start", "End", "GeneID")
    bids$Strand <- "+"
    bids[1:2,]
    write.table(bids[,c("GeneID", "Chr", "Start", "End", "Strand")], 
                paste0(out_dir, prefix, "_MUMERGE_tfit,dreg_", date, "_filt.saf"),
               quote=FALSE, sep="\t", row.names=FALSE)
    } 
if (type == "MU_COUNTS" | type == "BOTH") {
    # read in the closest bids overlapping file
    closest <- fread(close_bid_file)
    start.time <- Sys.time()
    # get the info needed for positions
    closest = get_info_for_calls(closest, count_win)
    # get the proper positions
    closest = get_new_pos(closest)
    # get the GTF lines and write them
    gtf_lines = get_GTF_lines(closest)
    writeLines(gtf_lines, paste0(out_dir, prefix, "_", date, "_mucounts_", str(count_win), ".gtf"))
    # now get the positive and negative forms
    gtf <- fread(paste0(out_dir, prefix, "_", date, "_mucounts_", str(count_win), ".gtf"))
    pos_gtf = gtf[gtf$V7 %in% c(".", "+"),]
    neg_gtf = gtf[gtf$V7 %in% c(".", "-"),]
    if (nrow(gtf)/3*2 != nrow(pos_gtf)) {
        stop("There was a problem in getting the GTF. One of the regions does not have transcripts on both strands.") }
    write.table(pos_gtf, paste0(out_dir, prefix, "_", date, "_mucounts_", str(count_win), "_pos.gtf"), 
            quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
    write.table(neg_gtf, paste0(out_dir, prefix, "_", date, "_mucounts_", str(count_win), "_neg.gtf"), 
            quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
    end.time <- Sys.time()
    cat("\nTime to get mu based counts:")
    print(end.time-start.time)
    }



##############################################
## 5. If TFEA output requested, print out ##
##############################################
if (TFEA == "YES") {
    # split up the tss and notss bids
    tss <- bids[bids$GeneID %in% over_filt$BidID,]
    nontss <- bids[!bids$GeneID %in% over_filt$BidID,]
    cat("\nNumber total Bids after filtering:", nrow(bids), "\n\tNumber TSS:", nrow(tss), "\tNumber NonTSS:", nrow(nontss))
    # save
    write.table(tss[,c("Chr", "Start", "End", "GeneID")], 
            paste0(out_dir, "tss_bid_", prefix, "_", date, "_forTFEA.bed"),
           quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
    write.table(nontss[,c("Chr", "Start", "End", "GeneID")], 
            paste0(out_dir, "nontss_bid_", prefix, "_", date, "_forTFEA.bed"),
           quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
    }
