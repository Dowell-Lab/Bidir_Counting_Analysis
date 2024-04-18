library(data.table)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
wd = args[1]
prefix = args[2]
date = args[3]
og_bid_file = args[4]
count_limit = args[5]

out_dir = paste0(wd, "regions/")
overlaps_dir = paste0(wd, "overlaps/")
tss_out_file = paste0(overlaps_dir, "overlaps_hg38_TSS1kb_withput_", prefix, "_MUMERGE_", date, ".bed")
count_out_file = paste0(overlaps_dir, "overlaps_hg38_withput_", prefix, "_MUMERGE_", date, ".bed")
closest_file = paste0(overlaps_dir, "closest_hg38_withput_", prefix, "_MUMERGE_", date, ".bed")

###################################
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
    # if more than one bidirectional has the same min distance --> hand annotate
    if (length(tss) != 1) {
        hand_annotate <- c(hand_annotate, isoform)
        cat("\nThe isoform", isoform, "has", length(tss), "TSS that are closest. Hand annotate")
    } else if (tss %in% TSS_used) {
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
count_over_neg = count_over_filt[count_over_filt$V6 == "-", ]
count_over_pos = count_over_filt[count_over_filt$V6 == "+", ]

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

##############################################
## 4. Make a SAF file of the remaining bids ##
##############################################
# read in the original bids
bids <- fread(og_bid_file)
pos_neg_remove <- union(pos_neg_colliding, pos_neg_filt)
cat("\nTotal Number of NonTSS bids removing due to convolution", length(pos_neg_remove))
bids <- bids[!bids$V4 %in% pos_neg_rem]
cat("\nTotal Number of Bids remaining", nrow(bids))
# save as a saf file
colnames(bids) <- c("Chr", "Start", "End", "GeneID")
bids$Strand <- "+"
bids[1:2,]
write.table(bids[,c("GeneID", "Chr", "Start", "End", "Strand")], 
            paste0(out_dir, prefix, "_MUMERGE_tfit,dreg_", date, "_filt.saf",
           quote=FALSE, sep="\t", row.names=FALSE)
