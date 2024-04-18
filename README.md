# Bidir_Counting_Analysis
This repo includes the generalizable pipeline for identification, counting, and length assessments of Bidirectionals from nascent RNA seq data


## bin
* 00_get_final_muMerge.sh
    * Inputs: mumerged separate Tfit and DReg bed files
    * Outputs: single file where Tfit (only those <2.5kb) and dREG regions are merged (if Tfit region overlaps dREG region, only keeps Tfit) and 4th column is the name in the format chr:start-stop-[tfit|dreg|tfit,dreg]
    * Analysis process: (uses get_final_muMerge.r)
        * First, remove Tfit regions >2.5kb since extremely low confidence
        * Second, only consider an overlap if >40% of either dREG or Tfit region is overlapping OR where the difference between mu calls <101bp (I've found this helps avoid merging calls that are actually different).
* 
