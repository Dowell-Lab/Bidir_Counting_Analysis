#!/bin/bash
#SBATCH --job-name=BidirCountFlow_Test
#SBATCH --output=/scratch/Users/hoto7260/nextflow_out/Bidir_Count/e_and_o/test_%j.out
#SBATCH --error=/scratch/Users/hoto7260/nextflow_out/Bidir_Count/e_and_o/test_%j.out
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --mem=1G
#SBATCH --partition short
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=hoto7260@colorado.edu

module load samtools/1.8
module load bedtools/2.28.0
module load openmpi/1.6.4
module load gcc/7.1.0
module load subread/1.6.2
module load R/4.4.0

source activate /Users/hoto7260/miniconda3/envs/python39
# python 3.9.17
# numpy  1.25.2
# pandas 2.0.3
 
 ## Clear the wokring directory
 work_dir=/scratch/Users/hoto7260/nextflow_out/Bidir_Count/tmp
 rm -rf ${work_dir}
 
 ## Run your code
NXF_DEBUG=1 /scratch/Shares/dowell/dbnascent/pipeline_assets/nextflow run /Users/hoto7260/src/Bidir_Counting_Analysis/main.nf \
 --mmfiltbams /scratch/Users/hoto7260/Bidir_Counting_Test/mmfiltbams/ \
 --cons_file /scratch/Users/hoto7260/Bidir_Counting_Test/cons_chr21_test.bed \
 -work-dir ${work_dir} \
 --outdir /scratch/Users/hoto7260/nextflow_out/Bidir_Count/Test/ \
 --prefix Testing \
 --date 2.10.25 \
 --count_type "Simple"
