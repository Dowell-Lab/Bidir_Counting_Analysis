#!/bin/bash

#SBATCH --job-name=Get_mmfilt_bams
#SBATCH --output=/scratch/Users/path/e_and_o/get_mmfiltbams_%j.out
#SBATCH --error=/scratch/Users/path/e_and_o/get_mmfiltbams_%j.out
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=50G
#SBATCH --partition short
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=email_address

crams=""
bams=""
# only need genome if starting with crams
genome=""

module load samtools

if crams != ""; do 
    for cram in ${crams}/*.sorted.cram; do
        prefix=$(basename "$cram" | sed 's/\.sorted\.cram$//')
        samtools view -@ 16 -b -1 -T ${genome} ${cram} > ${bams}/${prefix}.sorted.bam
        samtools index ${prefix}.sorted.bam ${prefix}.sorted.bam.bai

for bam in ${bams}/*.sorted.bam; do
    prefix=$(basename "$bam" | sed 's/\.sorted\.bam$//')
    samtools view -@ 32 -h -q 1 ${bam} | \
           grep -P '(NH:i:1|^@)' | \
           samtools view -h -b > ${bams}/${prefix}.mmfilt.sorted.bam
       samtools index ${bams}/${prefix}.mmfilt.sorted.bam ${bams}/${prefix}.mmfilt.sorted.bam.bai

done