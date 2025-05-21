library(data.table)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
wd = args[1]
prefix = args[2]
date = args[3]
count_req = args[4] # usually 10 for 4 samples
count_limit_genes = as.numeric(args[5])
count_type = args[6] # whetehr or not SIMPLE or MU_COUNTS
cat("\nARGUMENTS:", args, "\n")


region_dir = paste0(wd, "/regions/")
fixed_count_dir = paste0(wd, "/fixed_counts/")
overlaps_dir = paste0(wd, "/overlaps/")
count_out_file = paste0(overlaps_dir, "overlaps_hg38_withput_", prefix, "_MUMERGE_", date, ".bed")
dwn_out_file = paste0(overlaps_dir, "overlaps_hg38_withput_dwnstm_", prefix, "_MUMERGE_", date, ".bed")
closest_file = paste0(overlaps_dir, "closest_hg38_withput_", prefix, "_MUMERGE_", date, ".bed")
count_file = paste0(wd, "/counts/", prefix, "_str_put_genes.txt")
tss_file = paste0(region_dir, "tss_bid_", prefix, "_", date, ".txt")
count_dir <- paste0(wd, "/counts")
fixed_count_dir <- paste0(wd, "/fixed_counts")


if (count_type == "SIMPLE") {
    uns_counts_file = paste0(count_dir, "/", prefix, "_uns_bidirs.txt")
    pos_counts_file = paste0(count_dir, "/", prefix, "_pos_bidirs.txt")
    neg_counts_file = paste0(count_dir, "/", prefix, "_neg_bidirs.txt")
    full_counts_file = paste0(fixed_count_dir, "/", prefix, "_fixed_bids_", date, ".txt")
    } else {
    uns_counts_file = paste0(count_dir, "/", prefix, "_mucounts_str_bidirs.txt")
    pos_counts_file = paste0(count_dir, "/", prefix, "_mucounts_pos_bidirs.txt")
    neg_counts_file = paste0(count_dir, "/", prefix, "_mucounts_neg_bidirs.txt")
    full_counts_file = paste0(fixed_count_dir, "/", prefix, "_mucounts_fixed_bids_", date, ".txt")
    }




########################
## Read in the counts ##
########################
# read in counts
uns_counts <- fread(uns_counts_file)
pos_counts <- fread(pos_counts_file)
neg_counts <- fread(neg_counts_file)

#######################################################################
## Get the Bids that need to be fixed due to close/overlapping genes ##
#######################################################################
# read in the overlaps of Bids with transcripts (either full or downstream)
count_over <- fread(count_out_file)
dwn_over <- fread(dwn_out_file)

# Get the TSS bids
tss <- fread(tss_file)
tss <- unique(tss$BidID)

# read in the prel gene counts
counts <- fread(count_file)
colnames(counts) <- c("chrom", "start", "stop", "bidir", "dot", "strand", "count", "cov_bases", "length", "cov_frac")
# ensure that the Genes and overlaps have the same Geneid format
counts$Geneid <- str_split_fixed(counts$bidir, ",", 2)[,1]
# get the high counts
high_counts = counts[counts$cov_frac > count_limit_genes,]$Geneid
# only consider overlaps with genes with high counts
count_over <- count_over[count_over$V4 %in% high_counts]
dwn_over <- dwn_over[dwn_over$V4 %in% high_counts]

# Split the overlaps into positive and negative strands
count_over_neg = count_over[count_over$V6 == "-", ]
count_over_pos = count_over[count_over$V6 == "+", ]
dwn_over_neg = dwn_over[dwn_over$V6 == "-", ]
dwn_over_pos = dwn_over[dwn_over$V6 == "+", ]

# Get the NONTSS overlapping on one strand strands
bid_pos_conflict <- setdiff(union(count_over_pos$V10, dwn_over_pos$V10), tss)
bid_neg_conflict <- setdiff(union(count_over_neg$V10, dwn_over_neg$V10), tss)
remove_bids = intersect(bid_pos_conflict, bid_neg_conflict)

cat("\nTotal Number of NonTSS bids removing due to convolution", length(remove_bids))
cat("\nTotal regions considered", nrow(pos_counts), nrow(uns_counts), nrow(neg_counts))
pos_counts = pos_counts[!pos_counts$Geneid %in% remove_bids,]
uns_counts = uns_counts[!uns_counts$Geneid %in% remove_bids,]
neg_counts = neg_counts[!neg_counts$Geneid %in% remove_bids,]

cat("\nTotal regions considered", nrow(pos_counts), nrow(uns_counts), nrow(neg_counts))


##########################
## Get the Fixed Counts ##
##########################
# use a substring only found in the samples (not the other headings)
look_str <- "bam"

# get the pos counts of nontss bidirectionals w/ gene on neg strand
pos_counts <- data.frame(pos_counts[pos_counts$Geneid %in% bid_neg_conflict,])
pos_counts[1:2,]
colnames(pos_counts) <- colnames(uns_counts)
# multiply by 2
print(colnames(pos_counts)[grep(look_str, colnames(pos_counts))])
sample = pos_counts[,colnames(pos_counts)[grep(look_str, colnames(pos_counts))]]
sample = sample*2
sample[1:2,]
pos_counts <- cbind(pos_counts[,c("Geneid", "Chr", "Start", "End", "Strand", "Length")], sample)
pos_counts[1:2,]

# get the neg counts of bidirectionals w/ gene on pos strand
neg_counts <- data.frame(neg_counts[neg_counts$Geneid %in% bid_pos_conflict,])
colnames(neg_counts) <- colnames(uns_counts)
# multiply by 2
sample = neg_counts[,colnames(neg_counts)[grep(look_str, colnames(neg_counts))]]
sample = sample*2
neg_counts <- cbind(neg_counts[,c("Geneid", "Chr", "Start", "End", "Strand", "Length")], sample)

# Concatenate fixed counts to original counts of unoverlapping
# First have the unstranded counts not include the conflicting ones
uns_counts2 <- uns_counts[!uns_counts$Geneid %in% bid_neg_conflict,]
uns_counts2 <- uns_counts2[!uns_counts2$Geneid %in% bid_pos_conflict,]
uns_counts2[1:2,]

cat("\nSHAPES OF DATAFRAMES", dim(uns_counts2), dim(pos_counts), dim(neg_counts), "\n")

# Combine and save
full_counts <- rbind(uns_counts2, pos_counts, neg_counts)
dim(full_counts)
nrow(uns_counts2)+nrow(pos_counts)+nrow(neg_counts)
full_counts <- data.frame(full_counts)
colnames(full_counts) <- colnames(uns_counts)
full_counts[1:2,]
print(dim(full_counts))

# save full counts
write.table(full_counts, full_counts_file, 
            sep="\t", quote=FALSE, row.names=FALSE)

###################################################################
## Get the Bids that should be considered for fixed gene regions ##
###################################################################

### Only keep bidirectionals if they have more than X counts across all samples
full_nontss_counts <- full_counts[!full_counts$Geneid %in% tss,]
full_nontss_counts_moreX <- full_nontss_counts[rowSums(full_nontss_counts[,colnames(full_nontss_counts)[grep(look_str,colnames(full_nontss_counts))]])>count_req,]
dim(full_nontss_counts_moreX)
nrow(full_nontss_counts_moreX)/nrow(full_nontss_counts)

full_nontss_counts_moreX<- full_nontss_counts_moreX[,c("Chr", "Start", "End", "Geneid")]

write.table(full_nontss_counts_moreX, 
            paste0(region_dir, prefix, "_nontssbid_uniqueid_above", count_req, "_", date, ".bed"), 
            quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)


