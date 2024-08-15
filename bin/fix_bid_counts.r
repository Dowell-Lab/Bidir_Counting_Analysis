library(data.table)

args <- commandArgs(trailingOnly = TRUE)
wd = args[1]
prefix = args[2]
date = args[3]
count_req = args[4] # usually 10 for 4 samples
count_limit_genes = as.integer(args[5])
cat("\nARGUMENTS:", args, "\n")

region_dir = paste0(wd, "/regions/")
fixed_count_dir = paste0(wd, "/fixed_counts/")
overlaps_dir = paste0(wd, "/overlaps/")
count_out_file = paste0(overlaps_dir, "overlaps_hg38_withput_", prefix, "_MUMERGE_", date, ".bed")
closest_file = paste0(overlaps_dir, "closest_hg38_withput_", prefix, "_MUMERGE_", date, ".bed")
tss_file = paste0(region_dir, "tss_bid_", prefix, "_", date, ".txt")


###################################
## 1. Get the fixed counts ##
###################################
count_dir <- paste0(wd, "/counts")
fixed_count_dir <- paste0(wd, "/fixed_counts")

# read in counts
uns_counts <- fread(paste0(count_dir, "/", prefix, "_uns_bidirs.txt"))
pos_counts <- fread(paste0(count_dir, "/", prefix, "_pos_bidirs.txt"))
neg_counts <- fread(paste0(count_dir, "/", prefix, "_neg_bidirs.txt"))

#######################################################################
## Get the Bids that need to be fixed due to close/overlapping genes ##
#######################################################################
# get the genes with >40 counts
counts <- fread(paste0(count_dir, "/", prefix, "_str_put_genes.txt"))
count_matrix = data.frame(counts[,7:ncol(counts)])
rownames(count_matrix) <- counts$Geneid
count_matrix[1:2,]
high_counts = count_matrix[rowSums(count_matrix) > count_limit_genes,]
high_counts = rownames(high_counts)

# Get the TSS bids
tss <- fread(tss_file)
tss <- unique(tss$BidID)

# read in the overlaps with genes 
over <- fread(count_out_file)
# only keep those with significant counts
over <- over[over$V4 %in% high_counts]

# only keep closest  genes that are being transcribed significantly
closest <- fread(closest_file)
closest <- closest[closest$V4 %in% high_counts]
# only keep those that are within 5kb to the 3prime end 
# for positive genes, means V11 <5000 & positive 
closest <- closest[(closest$V6 == "+" & closest$V11 > 0 & closest$V11 < 5001) | closest$V6 == "-",]
# for negative genes, means V11 <5000 & negative
closest <- closest[(closest$V6 == "-" & closest$V11 < 0 & closest$V11 > -5001) | closest$V6 == "+",]

# only consider bids not previously filtered out
over <- over[over$V10 %in% uns_counts$Geneid]
closest <- closest[closest$V10 %in% uns_counts$Geneid]

# remove the tss and get the pos & neg ones
over <- over[!over$V10 %in% tss,]
over_bid_pos <- unique(over[over$V6 == "+",]$V10)
over_bid_neg <- unique(over[over$V6 == "-",]$V10)
closest <- closest[!closest$V10 %in% tss,]
closest_bid_pos <- unique(closest[closest$V6 == "+",]$V10)
closest_bid_neg <- unique(closest[closest$V6 == "-",]$V10)

# remove those that are within 5kb of 3' of transcribed gene & in opp strand gene
bid_pos_conflict <- union(closest_bid_pos, over_bid_pos)
bid_neg_conflict <- union(closest_bid_neg, over_bid_neg)
remove_bids <- intersect(bid_pos_conflict, bid_neg_conflict)
cat("\nBids to Remove", remove_bids[1:2])

# remove the bids that can't easily be deconvoluted
uns_counts <- uns_counts[!uns_counts$Geneid %in% remove_bids]
pos_counts <- pos_counts[!pos_counts$Geneid %in% remove_bids]
neg_counts <- neg_counts[!neg_counts$Geneid %in% remove_bids]

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

# save full counts
write.table(full_counts, 
            paste0(fixed_count_dir, "/", prefix, "_fixed_bids_", date, ".txt"), 
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


