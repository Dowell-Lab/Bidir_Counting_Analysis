## Hope Townsend 7/6/2023
## Script to get the regions of genes that are NOT overlapping a bidirectional & graph results

import numpy as np
import pandas as pd
import sys
import argparse
import func_get_fixed_gene_regions as fix
#import matplotlib.pylab as plt
import time

def main():
    
    ## GET ARGUMENTS
    # make parser to read arguments through command line
    parser = argparse.ArgumentParser()
    # add arguments
    parser.add_argument('--overlap_file', dest='overlaps', type=str)
    parser.add_argument('--use_bids_bedfile', dest='use_bids_file', type=str, default=None)
    parser.add_argument('--new_regions_output_file', dest='new_regions_output_file', type=str)
    parser.add_argument('--stats_output_file', dest='stats_output_file', type=str)
    #parser.add_argument('--stats_png_file', dest='stats_png_file', type=str)
    parser.add_argument('--full_gene_bed', dest='full_gene_bed', type=str, default="../data/processed_annotations/hg38_above850bplongisof_5ptrun_750dn_07_18_23.bed")
    parser.add_argument('--TSS_bid_file', dest='TSS_bid_file', type=str)
    # grab the arguments where arg is a dict
    arg = parser.parse_args()
    
    print("Arguments Used:", arg)
    
    ## CONVERT OVERLAP FILE TO APPROPRIATE FORM
    print("Reading in Overlaps")
    overlaps = pd.read_csv(arg.overlaps, sep="\t", header=None)
    # give it the appropriate column names
    if len(overlaps.columns) == 11:
        overlaps.columns = ['chr', 'Gene_Start', "Gene_Stop", "TranscriptID", "score", "strand", 
              "chr2", "Bid_Start", "Bid_Stop", "BidID", "overlap"]
    elif len(overlaps.columns) == 13:
        overlaps.columns = ['chr', 'Gene_Start', "Gene_Stop", "TranscriptID", "score", "strand", "GeneID", "Orig_Length",
              "chr2", "Bid_Start", "Bid_Stop", "BidID", "overlap"]
    else:
        print(overlaps.iloc[1:2,])
        raise ValueError("Not the proper number of columns")
    # FILTER OVERLAPS TO ONLY CONSIDER CERTAIN BIDIRECTIONALS IF DESIRED
    if len(arg.use_bids_file) is not None:
        # read in bed file to use
        old = overlaps.shape[0]
        use_bids = pd.read_csv(arg.use_bids_file, sep="\t", header=None)
        print("Only keeping", use_bids.shape[0], "bidirectionals to consider")
        use_bids.columns = ["chr", "start", "stop", "BidID"]
        overlaps = overlaps.loc[overlaps["BidID"].isin(use_bids.BidID)]
        print("Lost", old-overlaps.shape[0], "overlaps")
    
    # claculate the difference in overlap between length in gene
    overlaps["Length"] = overlaps.Gene_Stop - overlaps.Gene_Start + 1
    overlaps["DIFF"] = overlaps.Length - overlaps.overlap
    print(overlaps.head())

    ## GET NEW REGIONS AND STATS
    print("Getting Regions & Statistics")
    if arg.TSS_bid_file is not None:
        # read in TSS bids
        print("Getting TSS bidirectionals to not include in stats")
        TSS_bids = pd.read_csv(arg.TSS_bid_file, sep="\t")
        print("Number of TSS bids also found in overlaps file")
        TSS_bids = set(list(TSS_bids.BidID))
        print(len(TSS_bids), len(TSS_bids.intersection(set(overlaps["BidID"])))), 
        print("Getting nonoverlap regions and stats")
        t0=time.time()
        stats_df = fix.get_nonoverlap_GTF_genestats_focus(overlaps, list(TSS_bids))
        new_regions_df = None
        t1=time.time()
    else:
        print("No TSS bidirectionals specified to not include in stats")
        t0=time.time()
        [new_regions_df, stats_df] = fix.get_nonoverlap_GTF(overlaps)
        t1=time.time()
    print("Time to get nonoverlaps and/or stats:", t1-t0)
     # 1. pandas Dataframe in GTF where exons are the non overlapping Gene regions
        # 2. pandas Dataframe with the following columns:
            # 1. TranscriptID: Transcript (unique)
            # 2. Bidirectionals: All bidirectionals overlapping
            # 3. Num_Bids: # of bidirectionals overlapping
            # 4. BP_remain: # bp of gene that can still be counted
            # 5. Fraction_Overlap: Fraction of gene overlapping with bidirectionals
            # 6. Mean_Overlap: Mean of all overlaps
            # 7. Median_Overlap: Median of all overlaps

    ## SAVE NEW REGIONS & STATS
    if new_regions_df is not None:
        print("Adding non-overlapping genes to regions")
        ## ADD NONOVERLAPPING GENES
        # read in the full list of long genes
        full_genes = pd.read_csv(arg.full_gene_bed, 
                          sep="\t", header=None)
        if full_genes.shape[1] == 8:
            full_genes.columns = ["chr", "start", "stop", "TranscriptID", "score", "strand", "GeneID", "Length"]
        elif full_genes.shape[1] == 6:
            full_genes.columns = ["chr", "start", "stop", "TranscriptID", "score", "strand"]
        # get transcripts that don't have overlaps
        non_overlapping = set(full_genes.TranscriptID).difference(set(overlaps.TranscriptID))
        print("Number of Genes with no overlaps=", len(non_overlapping))
        # add these transcripts
        for transcript_name in non_overlapping:
            # get information
            filtered = full_genes[full_genes.TranscriptID == transcript_name]
            data = {"chr":list(filtered.chr)*2, "source":[".","."], 
                        "feature":["gene", "exon"],
                        "start": list(filtered.start)*2, 
                        "end": list(filtered.stop)*2,
                        "score":[".","."], "strand":list(filtered.strand)*2, "frame":[".","."], 
                        "attribute": ['gene_id '+transcript_name.split(":")[0], 'gene_id '+transcript_name.split(":")[0]+'; transcript_id '+transcript_name+'; region_id 1']}
            # add these data to the new regions gtf
            new_regions_df = pd.concat([new_regions_df,pd.DataFrame(data)])
    
        print("Shape of new regions file:", new_regions_df.shape)
        print("Saving New Regions")
        # ensure start and end are integers
        new_regions_df['start']=new_regions_df['start'].astype(int)
        new_regions_df['end']=new_regions_df['end'].astype(int)
        # ensure the proper numbers
        print("# of genes in new regions", new_regions_df[new_regions_df.feature == "gene"].shape[0])
        print("# of unique genes", len(set(full_genes.TranscriptID)))
        print(new_regions_df.head())
        new_regions_df.to_csv(arg.new_regions_output_file, sep="\t", index=False, header=False)
    print("Saving Stats")
    stats_df.to_csv(arg.stats_output_file, sep="\t", index=False)

#     ## GRAPH HISTOGRAMS of STATS
#     fig, axs = plt.subplots(ncols=1, nrows=5, figsize=(5.4, 8), 
#                             layout="constrained")
#     axs[0].hist(list(stats_df.Num_Bids), bins=50)
#     axs[0].set_title("# Bidirectionals per gene")
#     axs[1].hist(list(stats_df.BP_remain), bins=50)
#     axs[1].set_title("# BP remaining per gene")
#     axs[2].hist(list(stats_df.Fraction_Remain), bins=50)
#     axs[2].set_title("Fraction of Gene Remaining")
#     axs[3].hist(list(stats_df.BP_remain_100bp), bins=50)
#     axs[3].set_title("# BP remaining per gene after removing <100bp")
#     axs[4].hist(stats_df.Fraction_Remain_100bp, bins=50)
#     axs[4].set_title("Fraction of Gene Remaining after removing <100bp")
#     fig.savefig(arg.stats_png_file)
    
if __name__ == "__main__":
    main()


