# Bidir_Counting_Analysis
This repo includes the generalizable pipeline for identification, counting, and length assessments of Bidirectionals from nascent RNA seq data

## Variables/Parameters to choose:
* SRC = path to this repository (e.g. ~/Bidir_Counting_Analysis/)
* WD = directory where all output will be saved
* Count_Limit_Bids = the number required for a bidirectional to be considered impeding on gene transcription (I did 25 for 4 samples)
* Count_Limit_Genes = the number required for a gene to be considered transcribed enough to convolute bid counting (I did 40 for 4 samples)
* Prefix = string prefix to refer to project (e.g. Sasse2019_nascent)
* Date = string to refer to the date this pipeline is run

## Output (all in WD)
* regions/
    * Count regions for bidirectionals: ${prefix}_MUMERGE_tfit,dreg_${date}_filt.saf
    * TSS bidirectionals: tss_bid_${prefix}_${date}.txt
    * Count regions for genes: ${prefix}_new_regions_len_posneg_gene_filt_divconv_diff53prim_5ptrunc_${date}.gtf
* overlaps/
    * Overlaps with TSS 1kb regions (including putative refseq isoforms): overlaps_hg38_TSS1kb_withput_${prefix}_MUMERGE_${date}.bed
    * Overlaps with full isoforms (including putative refseq isoforms): overlaps_hg38_withput_${prefix}_MUMERGE_${date}.bed
    * Overlaps with the truncated refseq isoforms (which are used for counting): overlaps_hg38_trunc_${prefix}_MUMERGE_${date}.bed
    * Bedtools closest output for bidirectionals and genes: closest_hg38_withput_${prefix}_MUMERGE_${date}.bed
* counts/
    * Stranded counts for the putative genes: ${prefix}_str_put_genes.txt
    * Unstranded and stranded (both strands) counts for bidirectionals: ${prefix}_uns_bidirs.txt, ${prefix}_pos_bidirs.txt, ${prefix}_neg_bidirs.txt 
* fixed_counts/
    * fixed counts for the genes: ${prefix}_str_fixed_genes_${date}.txt
    * fixed counts for the bidirectionals: ${prefix}_MUMERGE_tfit,dreg_${date}_filt.saf 

## Pipeline Steps
* This pipeline starts with Tfit and dREG bed files that have been mumerged separately.
    * For Tfit, I recommend using the 3’ bedgraph (or equivalent) as Tfit is less likely to capture noise as being a bidirectional. Notes on how to determine the best bedgraph will eventually be added below. I also recommend that you ensure Tfit uses preliminary regions that include the TSS of genes AND the 800k regions called as having potential bidirectionals across multiple datasets. Details can be found below BUT this is automatically included in the Bidirectional Flow github repo on branch [Tfit_focus](https://github.com/Dowell-Lab/Bidirectional-Flow/tree/Tfit_focus). *NOTE*: The 3' methodology currently only works for single end (full warning in the README)
* 00_get_final_muMerge.sh
    * Inputs: mumerged separate Tfit and DReg bed files
    * Outputs: single file where Tfit (only those <2.5kb) and dREG regions are merged (if Tfit region overlaps dREG region, only keeps Tfit) and 4th column is the name in the format chr:start-stop-[tfit|dreg|tfit,dreg]
    * Analysis process: (uses get_final_muMerge.r)
        * First, remove Tfit regions >2.5kb since extremely low confidence
        * Second, only consider an overlap if >40% of either dREG or Tfit region is overlapping OR where the difference between mu calls <101bp (I've found this helps avoid merging calls that are actually different).
* 01a_Get_overlaps.sh
    * Inputs: Final bidirectionals used (from 00_get_final_muMerge.sh) with the windows for counting (I use 600bp) and TSS identification (I use 50bp)
    * Outputs: overlap files for future consideration of fixes AND to get the TSS bidirectionals
* 01b_Get_prel_gene_counts.sh
    * Inputs: Crams
    * Outputs: stranded counts of all gene isoforms
* 02_get_TSS_filt_bids.sh
    * Inputs: Nothing new (just direct to files made from previous scripts)
    * Outputs: SAF file for counting filtered Bids, TSS bids file
    * Important Note: I originally ran this with a jupyter notebook to allow more customizability. This version can be found under notebooks.
* 03_get_fix_bids_counts.sh
    * Inputs: SAF file from 02_get_TSS_filt_bids.sh, parameter options for # rowSum counts for gene to be considered disruptting Bid counts and vice versa
    * Outputs: Stranded & unstranded counts for filtered Bids (fixed and unfixed)
* 04_get_fix_gene_counts.sh
    * Inputs: Just edit names appropriately
    * Outputs: Fixed regions for Genes (GTF) and fixed counts for genes