#!/bin/bash

#SBATCH --job-name=count_gene
#SBATCH --output=/scratch/Users/hoto7260/Resp_Env/Sasse2019/e_and_o/featurecounts_gene1_%j.out
#SBATCH --error=/scratch/Users/hoto7260/Resp_Env/Sasse2019/e_and_o/featurecounts_gene1_%j.err
#SBATCH --time=05:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=5G
#SBATCH --partition short
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=hoto7260@colorado.edu

##load modules
module load samtools/1.8
module load subread/1.6.2

##########################
# EDIT THE FOLLOWING
##########################
## Naming
prefix="Sasse2019"
date=04_04_24
## OUT Directories
# wd is where the bams and counts will be saved
wd=/scratch/Users/hoto7260/Resp_Env/Sasse2019/out
bams=${wd}/bams
counts=${wd}/counts
fixed_counts=${wd}/fixed_counts
## IN Directories (CRAMS)
crams=/Shares/dbnascent/Sasse2019nascent/crams 

# genes
genes=/scratch/Shares/dowell/hoto7260/Bidir_Count/hg38_refseq_diff53prime_with_putatives_5ptrunc.saf

##########################
# CONVERT CRAMS TO BAMS
##########################

declare -a SRRs=(
[0]=SRR8429054
[1]=SRR8429055
[2]=SRR8429056
[3]=SRR8429057
[4]=SRR8429058
[5]=SRR8429059
)

for i in "${!SRRs[@]}"; do
    cram=${crams}/${SRRs[i]}.sorted.cram
    echo "Converting cram file to bam file for: " ${cram}
    samtools view -@ 8 -b -1 ${cram} > ${bams}/${SRRs[i]}.sorted.bam
    samtools index ${bams}/${SRRs[i]}.sorted.bam ${bams}/${SRRs[i]}.sorted.bam.bai
done

ls ${crams} | wc -l
ls ${bams} | wc -l


########################################## 
# Count reads over gene coordinates     ##
##########################################

echo "Submitted counts with Rsubread for original genes by strand......"
# Use this many threads (-T)
# Count multi-overlapping read
# Count reads in a strand specific manner (-s 2)
# Count by the gene feature (-t 'gene')
featureCounts \
    -T 32 \
    -O \
    -s 1 \
    -a ${genes} \
    -F 'SAF' \
    -o ${counts}/${prefix}_str_put_genes.txt \
    ${bams}/*.sorted.bam
    