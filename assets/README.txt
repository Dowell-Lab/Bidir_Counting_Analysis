Code to get the 10kb downstream regions (R 4.4):
    library(data.table)
    # get the full gene
    trunc <- fread("./assets/hg38_refseq_diff53prime_with_putatives_fixnames.sorted.bed")
    trunc[1:2,]
    # get the 10kb region downstream
    trunc$dwnstream_start = as.numeric(NA)
    trunc$dwnstream_end = as.numeric(NA)
    trunc[trunc$V6 == "+",]$dwnstream_start = trunc[trunc$V6 == "+",]$V3
    trunc[trunc$V6 == "-",]$dwnstream_end = trunc[trunc$V6 == "-",]$V2
    trunc[trunc$V6 == "+",]$dwnstream_end = trunc[trunc$V6 == "+",]$V3 + 10000
    trunc[trunc$V6 == "-",]$dwnstream_start = trunc[trunc$V6 == "-",]$V2 - 10000
    trunc[1:2,]
    dim(trunc[trunc$dwnstream_end < 1,]) # was 0
    write.table(trunc[,c("V1", "dwnstream_start", "dwnstream_end", "V4", "V5", "V6")], 
                "./assets/hg38_refseq_diff53prime_10kbdwnstm_with_putatives_fixnames.sorted.bed", 
               quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)