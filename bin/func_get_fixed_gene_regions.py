import pandas as pd
import numpy as np

def new_string(arg):
    return ":" + str(arg)

# function to get the new regions
def get_non_overlap_regions(gene_start, gene_stop, sorted_starts, sorted_ends):
    #INPUTS:
        # gene_start = int, starting coordinate of transcript body
        # gene_stop = int, ending coordinate of transcript body
        # sorted_starts = starts of bidirectionals sorted ascending
        # sorted_ends = ends of bidirectionals corresponding to same order as starts
    #OUTPUT:
        # new_region_list = list of tuples corresdponing to nonoverlapping regions
    i=0
    j=0
    gene_end = False
    new_region_list = []
    # == Get appropriate nonoverlapping start ==
    # if the bidirectional starts within gene body
    if sorted_starts[i] > gene_start:
        # the first region is at gene start
        new_region_list.append((gene_start, sorted_starts[i]))
    # == Get the appropriate end ==
    # if the bidirectional ends within gene body
    if sorted_ends[-1] < gene_stop:
        # then new regions should end with gene
        region_end = gene_stop
        gene_end = True
    # otherwise, the end should be the start of the last bidirectional
    else: 
        region_end = sorted_starts[-1]
    # == Get remaining regions ==
    # while the sorted starts < region_end
    for i in range(i, len(sorted_starts)-1):
        if sorted_starts[i+1] < region_end:
            new_region_list.append((sorted_ends[j], sorted_starts[i+1]))
            j=j+1
    # == Add the last region ==
    if gene_end: 
        j=-1
    else:
        j=-2
    new_region_list.append((sorted_ends[j], region_end))
    return new_region_list

def get_non_overlap_regions_gtf(transcript_name, chrom, strand, 
                               gene_start, gene_stop, sorted_starts, sorted_ends, 
                              diff_list, read_size=100):
    #INPUTS:
        # gene_start = int, starting coordinate of transcript body
        # gene_stop = int, ending coordinate of transcript body
        # sorted_starts = starts of bidirectionals sorted ascending
        # sorted_ends = ends of bidirectionals corresponding to same order as starts
        # diff_list = difference b/n length of gene & overlap of bid
    #OUTPUT:
        # new_region_df = dataframe with columns chr, source, feature, start, end, score, strand, frame, attribute
        #.     each of the regions where GeneID is exon 
    # IF ONE BID COMPLETELY OVERLAPS with less than read size available (75):
    if len(diff_list) ==1 and diff_list[-1] < read_size:
        data = {"chr":[str(chrom), str(chrom)], "source":[".","."], 
                "feature":["gene", "exon"],
                "start": [gene_start, gene_start], 
                "end": [gene_stop, gene_stop],
                "score":[".","."], "strand":[strand, strand], "frame":[".","."], 
                "attribute": ['gene_id '+str(transcript_name.split(":")[0]), 'gene_id '+str(transcript_name.split(":")[0])+'; transcript_id '+str(transcript_name) + '; region_id 1']}
        return pd.DataFrame(data)
    i=0
    j=0
    gene_end = False
    new_start_list = []
    new_end_list = []
    # == Get the appropriate end ==
    # if the last bidirectional ends within or at gene body
    if sorted_ends[-1] < gene_stop:
        # then new regions should end with gene
        region_end = gene_stop
        gene_end = True
    # otherwise, the end should be the start of the last bidirectional
    else: 
        region_end = sorted_starts[-1]
     
    # == Get appropriate nonoverlapping start ==
    # if the bidirectional starts AT the gene body
    if sorted_starts[i] == gene_start:
        # if only one region then use end of bidirectional and region_end
        if len(sorted_starts)==1:
            data = {"chr":[str(chrom), str(chrom)], "source":[".","."], 
                "feature":["gene", "exon"],
                "start": [gene_start, sorted_ends[-1]], 
                "end": [gene_stop, region_end],
                "score":[".","."], "strand":[strand, strand], "frame":[".","."], 
                "attribute": ['gene_id '+str(transcript_name.split(":")[0]), 'gene_id '+str(transcript_name.split(":")[0])+'; transcript_id '+str(transcript_name) + '; region_id 1']}
            return pd.DataFrame(data)
        
    # if the bidirectional starts within gene body
    elif sorted_starts[i] > gene_start:
        # the first region is at gene start
        new_start_list.append(gene_start)
        new_end_list.append(sorted_starts[i])
        # # if only one region 
        if len(sorted_starts)==1 and not gene_end:
            data = {"chr":[str(chrom), str(chrom)], "source":[".","."], 
                "feature":["gene", "exon"],
                "start": [gene_start, gene_start], 
                "end": [gene_stop, sorted_starts[i]],
                "score":[".","."], "strand":[strand, strand], "frame":[".","."], 
                "attribute": ['gene_id '+str(transcript_name.split(":")[0]), 'gene_id '+str(transcript_name.split(":")[0])+'; transcript_id '+str(transcript_name) + '; region_id 1']}
            return pd.DataFrame(data)
    # if only one bidirectional that starts before gene body
    elif len(sorted_starts)==1:
        # get end and gene start
        data = {"chr":[str(chrom), str(chrom)], "source":[".","."], 
                "feature":["gene", "exon"],
                "start": [gene_start, sorted_ends[-1]], 
                "end": [gene_stop, region_end],
                "score":[".","."], "strand":[strand, strand], "frame":[".","."], 
                "attribute": ['gene_id '+str(transcript_name.split(":")[0]), 'gene_id '+str(transcript_name.split(":")[0])+'; transcript_id '+str(transcript_name) + '; region_id 1']}
        return pd.DataFrame(data)
        # == Get remaining regions ==
    # while the sorted starts < region_end
    for i in range(i, len(sorted_starts)-1):
        if sorted_starts[i+1] < region_end:
            new_start_list.append(sorted_ends[j])
            new_end_list.append(sorted_starts[i+1])
            j=j+1
    # == Add the last region ==
    # if the gene is the end then use last bidirectional's end, otherwise 2nd to last
    if gene_end: 
        j=-1
    else:
        j=-2
    new_start_list.append(sorted_ends[j])
    new_end_list.append(region_end)
    # get the # of regions
    num = len(new_start_list)
    def gtf_string(arg):
        return str('gene_id ' + str(transcript_name.split(":")[0])+'; transcript_id ' + str(transcript_name) + '; region_id ' + str(arg))
    GeneID = list(map(gtf_string, range(1,num+1)))
    data = {"chr":[str(chrom)]*(num+1), "source":["."]*(num+1), 
                "feature":["gene"]+["exon"]*num,
                "start": [gene_start]+new_start_list, 
                "end": [gene_stop]+new_end_list,
                "score":["."]*(num+1), "strand":[strand]*(num+1), "frame":["."]*(num+1), 
                "attribute": ['gene_id '+str(transcript_name.split(":")[0])] + GeneID}
    return pd.DataFrame(data)

# function to get the new regions
def get_non_overlap_regions_df(transcript_name, chrom, strand, 
                               gene_start, gene_stop, sorted_starts, sorted_ends, 
                              diff_list):
    #INPUTS:
        # gene_start = int, starting coordinate of transcript body
        # gene_stop = int, ending coordinate of transcript body
        # sorted_starts = starts of bidirectionals sorted ascending
        # sorted_ends = ends of bidirectionals corresponding to same order as starts
        # diff_list = difference b/n length of gene & overlap of bid
    #OUTPUT:
        # new_region_df = dataframe with columns GeneID, Chr, Start, End, Strand for 
        #.     each of the regions where GeneID is transcript+number
    # IF BID COMPLETELY OVERLAPS with only 10bp remaining:
    if diff_list[-1] < 10:
        data = {"GeneID": [transcript_name + ":BID"], 
                "Chr": [str(chrom)], 
                "Start": gene_start, 
                "End": gene_stop, 
                "Strand": [strand]}
        return pd.DataFrame(data)
    i=0
    j=0
    gene_end = False
    new_start_list = []
    new_end_list = []
    # == Get the appropriate end ==
    # if the last bidirectional ends within or at gene body
    if sorted_ends[-1] < gene_stop:
        # then new regions should end with gene
        region_end = gene_stop
        gene_end = True
    # otherwise, the end should be the start of the last bidirectional
    else: 
        region_end = sorted_starts[-1]
     
    # == Get appropriate nonoverlapping start ==
    # if the bidirectional starts AT the gene body
    if sorted_starts[i] == gene_start:
        # if only one region then use end of bidirectional and region_end
        if len(sorted_starts)==1:
            data = {"GeneID": [transcript_name + ":1"], 
                "Chr": [str(chrom)], 
                "Start": sorted_ends[-1], 
                "End": region_end, 
                "Strand": [strand]}
            return pd.DataFrame(data)
        
    # if the bidirectional starts within gene body
    elif sorted_starts[i] > gene_start:
        # the first region is at gene start
        new_start_list.append(gene_start)
        new_end_list.append(sorted_starts[i])
        # # if only one region 
        if len(sorted_starts)==1 and not gene_end:
            data = {"GeneID": [transcript_name + ":1"], 
                "Chr": [str(chrom)], 
                "Start": new_start_list, 
                "End": new_end_list, 
                "Strand": [strand]}
            return pd.DataFrame(data)
    # if only one bidirectional that starts before gene body
    elif len(sorted_starts)==1:
        # get end and gene start
        data = {"GeneID": [transcript_name + ":1"], 
                "Chr": [str(chrom)], 
                "Start": sorted_ends[-1],
                "End": region_end, 
                "Strand": [strand]} # originally start=gene_start
        return pd.DataFrame(data)
   
        
    # == Get remaining regions ==
    # while the sorted starts < region_end
    for i in range(i, len(sorted_starts)-1):
        if sorted_starts[i+1] < region_end:
            new_start_list.append(sorted_ends[j])
            new_end_list.append(sorted_starts[i+1])
            j=j+1
    # == Add the last region ==
    # if the gene is the end then use last bidirectional's end, otherwise 2nd to last
    if gene_end: 
        j=-1
    else:
        j=-2
    new_start_list.append(sorted_ends[j])
    new_end_list.append(region_end)
    # get the # of regions
    num = len(new_start_list)
    # get gene ids
    add = list(map(new_string, range(1,num+1)))
    GeneID = list(map(lambda x: transcript_name + x, add))
    data = {"GeneID": GeneID, 
            "Chr": [str(chrom)]*num, 
            "Start": new_start_list, 
            "End": new_end_list, 
            "Strand": list(strand)*num}
    return pd.DataFrame(data)
    

def get_nonoverlap_SAF(overlaps, TSS_bids=None):
    # INPUTS: overlaps: dataframe including columns TranscriptID, chr, Length, Gene_Start, Gene_Stop, strand, BidID, Bid_Start, Bid_Stop, overlap
    # OUTPUTS: list of
        # 1. pandas Dataframe in SAF format of the non overlapping Gene regions
        # 2. pandas Dataframe with the following columns:
            # 1. TranscriptID: Transcript (unique)
            # 2. Bidirectionals: All bidirectionals overlapping
            # 3. Num_Bids: # of bidirectionals overlapping
            # 4. BP_remain: # bp of gene that can still be counted
            # 5. Fraction_Overlap: Fraction of gene overlapping with bidirectionals
            # 6. Mean_Overlap: Mean of all overlaps
            # 7. Median_Overlap: Median of all overlaps
    transcript_list = list(set(overlaps["TranscriptID"]))
    bid_list = [] # bidirectionals
    num_bids = [] # num bidirectionals overlapping
    per_overlap_list = [] # percentage overlap
    num_bp_remain_list = [] # num_bp remaining
    mean_overlap_list = [] # average overlap
    median_overlap_list = [] # median overlap
    new_regions_df = pd.DataFrame({"GeneID":[], "Chr":[], 
                                   "Start":[], "End":[], "Strand":[]}) # list of regions
    # if not supposed to include TSSs
    if TSS_bids is not None:
        bid_list_notss = []
        num_bids_notss = []
    # for each transcript ID
    for transcript in transcript_list:
        print("TRANSCRIPT", transcript)
        # Get the bidirectionals overlapping
        filtered = overlaps.loc[overlaps.TranscriptID == transcript,]
        # sort based on start
        filtered = filtered.sort_values(by=['Bid_Start'])
        # get the bidirectionals & #
        bid_list.append(list(filtered.BidID))
        num_bids.append(len(filtered.BidID))
        # if also supposed to establish for non TSS
        if TSS_bids is not None:
            # see which bidirectionals NOT in TSS. 
            notss_bids = set(filtered.BidID).difference(TSS_bids)
            bid_list_notss.append(list(notss_bids))
            num_bids_notss.append(len(notss_bids))
        # get the mean & median overlap
        mean_overlap_list.append(np.mean(filtered.overlap))
        median_overlap_list.append(np.median(filtered.overlap))
        # get the # bp remaining
        bp_taken = sum(list(filtered.overlap))
        num_bp_remain_list.append(list(filtered.Length)[0] - bp_taken)
        per_overlap_list.append(round(bp_taken/ list(filtered.Length)[0], 3))
        # get the new regions
        df = get_non_overlap_regions_df(transcript, list(filtered.chr)[0], list(filtered.strand)[0], 
                                        list(filtered.Gene_Start)[0], list(filtered.Gene_Stop)[0], 
                                        list(filtered.Bid_Start), list(filtered.Bid_Stop), list(filtered.DIFF))
        new_regions_df = pd.concat([new_regions_df,df])
    if TSS_bids is None:
        stat_df = {"TranscriptID":transcript_list,  "Bidirectionals":bid_list, 
                   "Num_Bids":num_bids, "BP_remain":num_bp_remain_list, 
                   "Fraction_Overlap":per_overlap_list, 
                   "Mean_Overlap":mean_overlap_list, "Median_Overlap":median_overlap_list}
    else:
        stat_df = {"TranscriptID":transcript_list,  "Bidirectionals":bid_list, 
                   "Num_Bids":num_bids, "BP_remain":num_bp_remain_list, 
                   "Fraction_Overlap":per_overlap_list, 
                   "Mean_Overlap":mean_overlap_list, "Median_Overlap":median_overlap_list, 
                  "Bidirectionals_nonTSS":bid_list_notss, "Num_Bids_nonTSS":num_bids_notss}
    stat_df = pd.DataFrame(stat_df) 
    return [new_regions_df, stat_df]

def get_nonoverlap_GTF(overlaps, TSS_bids=None):
    # INPUTS: overlaps: dataframe including columns TranscriptID, chr, Length, Gene_Start, Gene_Stop, strand, BidID, Bid_Start, Bid_Stop, overlap
    # OUTPUTS: list of
        # 1. pandas Dataframe in GTF format of the non overlapping Gene regions:
            # chr, source, feature, start, end, score, strand, frame, attribute where each region is consideered an
            #     exon of the gene
        # 2. pandas Dataframe with the following columns:
            # 1. TranscriptID: Transcript (unique)
            # 2. Bidirectionals: All bidirectionals overlapping
            # 3. Num_Bids: # of bidirectionals overlapping
            # 4. BP_remain: # bp of gene that can still be counted
            # 5. Fraction_Overlap: Fraction of gene overlapping with bidirectionals
            # 6. Mean_Overlap: Mean of all overlaps
            # 7. Median_Overlap: Median of all overlaps
    transcript_list = list(set(overlaps["TranscriptID"]))
    bid_list = [] # bidirectionals
    num_bids = [] # num bidirectionals overlapping
    per_remain_list = [] # percentage remaining
    num_bp_remain_list = [] # num_bp remaining
    num_bp_remain_post_list = [] # num_bp remaining after removing <100bp regions
    per_remain_post_list = [] # percentage remaining after removing <100bp regions
    new_regions_df = pd.DataFrame({"chr":[], "source":[], 
                "feature":[],"start": [], "end": [],
                "score":[], "strand":[], "frame":[], 
                "attribute": []}) # list of regions
    # if not supposed to include TSSs
    if TSS_bids is not None:
        bid_list_notss = []
        num_bids_notss = []
    # for each transcript ID
    for transcript in transcript_list:
        #print("TRANSCRIPT", transcript)
        # Get the bidirectionals overlapping
        filtered = overlaps.loc[overlaps.TranscriptID == transcript,]
        # sort based on start
        filtered = filtered.sort_values(by=['Bid_Start'])
        # print("FILLTERED")
        # print(filtered)
        # get the bidirectionals & #
        bid_list.append(list(filtered.BidID))
        num_bids.append(len(filtered.BidID))
        # if also supposed to establish for non TSS
        if TSS_bids is not None:
            # see which bidirectionals NOT in TSS. 
            notss_bids = set(filtered.BidID).difference(TSS_bids)
            bid_list_notss.append(list(notss_bids))
            num_bids_notss.append(len(notss_bids))
        
        # get the new regions
        df = get_non_overlap_regions_gtf(transcript, list(filtered.chr)[0], list(filtered.strand)[0], 
                                        list(filtered.Gene_Start)[0], list(filtered.Gene_Stop)[0], 
                                        list(filtered.Bid_Start), list(filtered.Bid_Stop), list(filtered.DIFF))
        # see if any regions below 100bp
        df['Length'] = df.end-df.start+1
        # print("DF")
        # print(df)
        # get the # bp remaining
        bp_remain = sum(list(df.Length))-list(filtered.Length)[0]
        num_bp_remain_list.append(bp_remain)
        per_remain_list.append(round((bp_remain/list(filtered.Length)[0]),3))
        df = df[df.Length >100]
        # add the new remaining bp and percentage
        bp_remain = sum(list(df.Length))-list(filtered.Length)[0]
        num_bp_remain_post_list.append(bp_remain)
        per_remain_post_list.append(round((bp_remain/list(filtered.Length)[0]),3))
        new_regions_df = pd.concat([new_regions_df,df.loc[:,df.columns!='Length']])
        # print("neew regions")
        # print(new_regions_df)
    if TSS_bids is None:
        stat_df = {"TranscriptID":transcript_list,  "Bidirectionals":bid_list, 
                   "Num_Bids":num_bids, "BP_remain":num_bp_remain_list, 
                   "Fraction_Remain":per_remain_list, 
                   "BP_remain_100bp":num_bp_remain_post_list, "Fraction_Remain_100bp":per_remain_post_list}
    else:
        stat_df = {"TranscriptID":transcript_list,  "Bidirectionals":bid_list, 
                   "Num_Bids":num_bids, "BP_remain":num_bp_remain_list, 
                   "Fraction_Remain":per_remain_list, 
                   "BP_remain_100bp":num_bp_remain_post_list, "Fraction_Remain_100bp":per_remain_post_list, 
                   "Bidirectionals_nonTSS":bid_list_notss, "Num_Bids_nonTSS":num_bids_notss}
    stat_df = pd.DataFrame(stat_df)
    print("Min Fraction remain", max(stat_df.Fraction_Remain))
    print("Min BP_remain", min(stat_df.BP_remain))
    print("Min Fraction main post removing regions <100bp", max(stat_df.Fraction_Remain_100bp))
    print("Min BP_remain post removing regions <100bp", min(stat_df.BP_remain_100bp))
    return [new_regions_df, stat_df]


def get_nonoverlap_GTF_genestats_focus(overlaps, TSS_bids=None):
    ## Same idea as get_nonoverlap_GTF but streamlined for full genes (so don't care about Fraction_Remain_100bp and no new regions file)
    # INPUTS: overlaps: dataframe including columns TranscriptID, chr, Length, Gene_Start, Gene_Stop, strand, BidID, Bid_Start, Bid_Stop, overlap
    # OUTPUTS: 
        # pandas Dataframe with the following columns:
            # 1. TranscriptID: Transcript (unique)
            # 2. Bidirectionals: All bidirectionals overlapping
            # 3. Num_Bids: # of bidirectionals overlapping
            # 4. BP_remain: # bp of gene that can still be counted
            # 5. Fraction_Remain: Fraction of gene not overlapped with bidirectionals
            # 6. Bidirectionals_nonTSS: Bidirectionals overlapping not including TSS one
            # 7. Num_Bids_nonTSS: # of bidirectionals overlapping not including TSS one
            # 8. BP_remain_nonTSS: # bp of gene without non-TSS bidirectionals overlapped
            # 9. Fraction_Remain_nonTSS: Fraction of gene without non-TSS bidirectionals overlapped
    transcript_list = list(set(overlaps["TranscriptID"]))
    print("Number of transcripts", len(transcript_list))
    bid_list = [] # bidirectionals
    num_bids = [] # num bidirectionals overlapping
    per_remain_list = [] # percentage remaining
    num_bp_remain_list = [] # num_bp remaining
    per_remain_notss_list = [] # percentage remaining if not counting TSS bidirectional
    num_bp_remain_notss_list = [] # num_bp remeaining if not counting TSS bidirectional
    bid_list_notss = []
    num_bids_notss = []
    # for each transcript ID
    for transcript in transcript_list:
        #print("TRANSCRIPT", transcript)
        # Get the bidirectionals overlapping
        filtered = overlaps.loc[overlaps.TranscriptID == transcript,]
        # sort based on start
        filtered = filtered.sort_values(by=['Bid_Start'])
        # print("FILLTERED")
        # print(filtered)
        # get the bidirectionals & #
        bid_list.append(list(filtered.BidID))
        num_bids.append(len(filtered.BidID))
        # if also supposed to establish for non TSS
        # see which bidirectionals NOT in TSS. 
        notss_bids = set(filtered.BidID).difference(TSS_bids)
        bid_list_notss.append(list(notss_bids))
        num_bids_notss.append(len(notss_bids))
        ## NUMBER BP remaining
        gene_length = list(filtered.Length)[0]
        # Regardless of TSS
        bp_remain = gene_length - sum(list(filtered.overlap))# Gene Length - overlaps (only works cuz no overlapping bids)
        num_bp_remain_list.append(bp_remain)
        per_remain_list.append(round(bp_remain/gene_length,3))
        # Considering TSS to not include in stats
        #   only get overlaps that don't include TSS bidirectionals
        filtered = filtered[~filtered["BidID"].isin(TSS_bids)]
        bp_remain = gene_length - sum(list(filtered.overlap))# Gene Length - overlaps (only works cuz no overlapping bids)
        num_bp_remain_notss_list.append(bp_remain)
        per_remain_notss_list.append(round(bp_remain/gene_length,3))
    print(len(transcript_list), len(bid_list), len(num_bids), len(num_bp_remain_list), 
          len(per_remain_list), len(per_remain_list), len(bid_list_notss), len(num_bids_notss), 
          len(num_bp_remain_notss_list), len(per_remain_notss_list))
    stat_df = {"TranscriptID":transcript_list,  "Bidirectionals":bid_list, 
                   "Num_Bids":num_bids, "BP_remain":num_bp_remain_list, 
                   "Fraction_Remain":per_remain_list,
                   "Bidirectionals_nonTSS":bid_list_notss, "Num_Bids_nonTSS":num_bids_notss, 
              "BP_remain_nonTSS":num_bp_remain_notss_list, "Fraction_Remain_nonTSS":per_remain_notss_list}
    stat_df = pd.DataFrame(stat_df)
    print("Min Fraction remain", max(stat_df.Fraction_Remain))
    print("Min BP_remain", min(stat_df.BP_remain))
    print("Min Fraction remain no TSS", max(stat_df.Fraction_Remain_nonTSS))
    print("Min BP_remain no TSS", min(stat_df.BP_remain_nonTSS))
    return stat_df

def get_overlap_bid_stats(overlaps):
    # SUMMARY: THis function looks at the overlaps from a bidirectional-centered view to get stats as described in "OUTPUTS."
    # INPUTS: overlaps: dataframe including columns TranscriptID, chr, Length, Gene_Start, Gene_Stop, strand, BidID, Bid_Start, Bid_Stop, Bid_Length, Frac_Over_Gene, Frac_Over_Bid, overlap
    # OUTPUTS: pandas Dataframe with the following columns:
            # 1. BidID: unique ID for bidirectional
            # 2. Transcripts: All gene_transcripts overlapping bidirectional (comma separated)
            # 3. Num_Transcripts: # of transcripts bidirectional overlaps
            # 4. Fraction_Overlaps_bid: Fractions of bidirectional overlapping with genes (comma separated)
            # 5. Fraction_Overlaps_genes: Fractions of genes overlapping with bidirectional (comma separated)
            # 6. Tot_Fraction_Overlap_bid: Total Fraction of bidirectional overlapping with genes
            
    bid_list = list(set(overlaps["BidID"])) # list of the unique bidirectionals
    trans_list = [] # transcripts overlapping
    num_trans = [] # Num transcripts overlapping
    strand_trans = [] # strands of transcripts overlapping
    frac_overlap_bid = [] # fractions of bid overlapping w/ genes
    frac_overlap_gene = [] # fractions of gene overlapping w/ bid
    frac_overlap_list = [] # total fraction of overlap of bid w/ genes

    # for each transcript ID
    for bid in bid_list:
        print("Bidirectional", bid)
        # Get the transcripts overlapping
        filtered = overlaps.loc[overlaps.BidID == bid,]
        # get the transcripts & num
        trans_list.append(list(filtered.TranscriptID))
        num_trans.append(len(filtered.TranscriptID))
        strand_trans.append(list(filtered.strand))
        # get the fraction overlaps
        frac_overlap_bid.append(list(filtered.Frac_Over_Bid))
        frac_overlap_gene.append(list(filtered.Frac_Over_Gene))
        # get the fraction of bidirectional overlappinig
        frac_overlap_list.append(sum(list(filtered.Frac_Over_Bid)))
    
    # get final statistics
    stat_df = {"BidID":bid_list,  "Transcripts":trans_list, 
                   "Num_Transcripts":num_trans, "Strand_Trans":strand_trans,
               "Fraction_Overlaps_bid":frac_overlap_bid, 
                   "Fraction_Overlaps_genes":frac_overlap_gene, 
                   "Tot_Fraction_Overlap_bid":frac_overlap_list}
    stat_df = pd.DataFrame(stat_df) 
    return stat_df
    