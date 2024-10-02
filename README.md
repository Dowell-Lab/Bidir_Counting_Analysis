# Bidir_Counting_Analysis
This repo includes the generalizable pipeline for identification of bidirectionals (TSS & NonTSS), and counting (while addressing overlapping transcription) with nascent RNA seq data

## Requirements:
* bedtools (2.28.0 used)
* samtools (1.8 used)
* R (3.6 - 4.4 tested)
    * data.table
    * stringR
    * If on a supercomputer (e.g. Fiji), you can get data.table by running the code below. You'll get an alert that the opt/R is not writable so it must install data.table into a personal R directory that it will store in a location it provides (usually ~/R/x86_64-pc-linux-gnu-library or some variation). You can say YES to this. From now on, when you install packages, they will be stored here. NOW, you can install stringR.
 ```
module load R/4.4.0
R
>install.packages("data.table")
```
* python >3 (3.6 & 3.9 tested)
    * pandas (2.2.2 tested with python 3.9)
    * numpy (1.26.5 tested with python 3.9)
 
-------------------------------------------------------------

## How to Run
* combined_flow.sh will run all steps of the pipeline. It will get the TSS/NonTSS bidirectionals from a file of consensus bidirectionals (usually from muMerge). Then, it will get the fixed counts of both bidirectionals and genes according to overlapping transcription. Details of the pipeline are in a later section. The output is below.
* Individual scripts for each step can be found in bin/
* If you just want to get TSS/NonTSS bidirectionals, skip to that appropriate section below.

### Variables/Parameters to choose:
* SRC = path to this repository (e.g. ~/Bidir_Counting_Analysis/)
* WD = directory where all output will be saved
* COUNT_WIN = This is the window used to make the bid counting file (new_start=mu-COUNT_WIN & new_stop=mu+COUNT_WIN).
    * Recommended = 400-700bp --> 800bp-1.4kb total regions
    * If using muMerge to make the consensus regions, the output lengths are a confidence interval around mu and can range a lot (e.g. 20bp-3kb). With low counts of bidirectionals, we want to ensure we are capturing enough of the region for maximum info while also avoiding improper noise. Therefore, we transform the muMerge lengths into fixed ones around mu.
* TSS_WIN = This is the window used to make the bid file for identifying TSS bidirectionals
    * Works the same as COUNT_WIN bu
* Count_Limit_Bids = the number required for a bidirectional to be considered impeding on gene transcription (I did 25 for 4 samples)
* Count_Limit_Genes = the number required for a gene to be considered transcribed enough to convolute bid counting (I did 40 for 4 samples)
* Prefix = string prefix to refer to project (e.g. Sasse2019_nascent)
* Date = string to refer to the date this pipeline is run
* COUNT_TYPE = string to refer to the type of counting method you want (options are "SIMPLE", "MU_COUNTS", or "BOTH" and are described below):
     * The "SIMPLE" approach uses a fixed window (COUNT_WIN) around the centers of the consensus regions providing without worrying about neighboring bidirectionals, and considering the entire region for positive and negative counts (doesn't split into transcripts).
     * The "MU_COUNTS" approach uses a fixed window (COUNT_WIN) for the region of the tRE BUT has two additions (Image below to help visualize):
          1) Considers overlapping regions: If a fixed counting window means the given regions are now overlapping, the neighboring tREs will be counted so that the maximum distance of a transcript is to the middle (mu) of the nearest neighboring tRE.
          2) Considers transcripts of each region separately: Assuming each consensus region is a bidirectional (transcripts from both strands originating at the middle of the region (mu), we then count each transcript separate from one another. This allows more exact deconvolution of overlapping bidirectionals while maintaining as much data as possible.
      * The "BOTH" option provides the counts and regions files for both of these methods.
* TFEA = string (YES or NO) of whether or not you want TSS & NonTSS bed files to be used as input into TFEA
* MAKE_NAMES = string (YES or NO) of whether or not you want the code to add unique names for the bidirectionals in CONS_BID. If NO, it will assume the fourth column is the unique names. Unique names are in the format chrom:start-stop.

### Just get the TSS/NonTSS bidirectionals to run TFEA
* get_TSS_nonTSSbids_forTFEA.sh will only run the first few steps of combined_flow.sh to get the TSS & NonTSS bidirectionals as input for TFEA or for counting. It will NOT do fixed counting for the regions but WILL provide counts (so no fixed_counts output).

-------------------------------------------------------------

## Input & Output
### Input Files Required
* Consensus bidirectionals (CONS_BID)
    * This should be the consensus regions of the bidirectionals in bed format (example is muMerge output). 
### Output (all in the working directory parameter chosen)
* regions/
    * Count window bed file: `${prefix}_${date}_MUMERGE_${count_win}win_TSS.sorted.bed`
    * Count regions for bidirectionals without convolution on both strands:
        * `${prefix}_MUMERGE_tfit,dreg_${date}_filt.saf` (**if COUNT_TYPE==SIMPLE or BOTH**)
        * `${prefix}_mucounts_${count_win}_${date}.gtf` (both transcripts), `${prefix}_mucounts_${count_win}_${date}_pos.gtf` (just positive transcripts), `${prefix}_mucounts_${count_win}_${date}_neg.gtf` (just negative transcripts) (**if COUNT_TYPE==MU_COUNTS or BOTH**)
    * TSS bidirectionals: `tss_bid_${prefix}_${date}.txt`
    * Count regions for genes: `${prefix}_new_regions_len_posneg_gene_filt_divconv_diff53prim_5ptrunc_${date}.gtf`
    * TSS bidirectionals bed (for TFEA): tss_bid_${prefix}_${date}_forTFEA.bed (**if TFEA**)
    * NonTSS bidirectionals bed (for TFEA): nontss_bid_${prefix}_${date}_forTFEA.bed (**if TFEA**)
* overlaps/
    * Overlaps with TSS 1kb regions (including putative refseq isoforms): `overlaps_hg38_TSS1kb_withput_${prefix}_MUMERGE_${date}.bed`
    * Overlaps with full isoforms (including putative refseq isoforms): `overlaps_hg38_withput_${prefix}_MUMERGE_${date}.bed`
    * Overlaps with the truncated refseq isoforms (which are used for counting): `overlaps_hg38_trunc_${prefix}_MUMERGE_${date}.bed`
    * Bedtools closest output for bidirectionals and genes: `closest_hg38_withput_${prefix}_MUMERGE_${date}.bed`
* counts/
    * Stranded counts for the putative genes: `${prefix}_str_put_genes.txt`
    * Unstranded and stranded (both strands) counts for bidirectionals: `${prefix}_uns_bidirs.txt`, `${prefix}_pos_bidirs.txt`, `${prefix}_neg_bidirs.txt`
* fixed_counts/
    * fixed counts for the genes: `${prefix}_str_fixed_genes_${date}.txt`
    * fixed counts for the bidirectionals: `${prefix}_MUMERGE_tfit,dreg_${date}_filt.saf` 

-------------------------------------------------------------

## Full Pipeline Details
* This pipeline starts with Tfit and dREG bed files that have been mumerged separately.
    * For Tfit, I recommend using the 3’ bedgraph (or equivalent) as Tfit is less likely to capture noise as being a bidirectional. Notes on how to determine the best bedgraph will eventually be added below. I also recommend that you ensure Tfit uses preliminary regions that include the TSS of genes AND the 800k regions called as having potential bidirectionals across multiple datasets. Details can be found below BUT this is automatically included in the Bidirectional Flow github repo on branch [Tfit_focus](https://github.com/Dowell-Lab/Bidirectional-Flow/tree/Tfit_focus). *NOTE*: The 3' methodology currently only works for single end (full warning in the README)
* To run the full pipeline with one script, you can use the combined_flow.sh
    * WARNINGS: 
        * I ran /bin/00_get_final_muMerge.sh separately when creating this file so feel free to run that separately or add the code
        * It is currently designed to work with bam files. You can refer to 01b_Get_prel_gene_counts.sh to add the code to convert the crams to bams. I can add this later.
* 00_get_final_muMerge.sh
    * Inputs: mumerged separate Tfit and DReg bed files, and option if want to keep sample id information from muMerge 
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
