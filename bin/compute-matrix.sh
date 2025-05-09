#!/bin/bash

#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10gb
#SBATCH --job-name=matrix
#SBATCH --error=./outputs/logs/matrix.err
#SBATCH --error=./outputs/logs/matrix_log.out

# This script uses a bed file with regions of interest and sample signal 
# intensity bigWig files to generate a signal intensity matrix for each
# sample. The matrix will show signal intensity in a sample spaning 1000
# base pairs before and after the center of each region of interest, in
# 10 base pair intervals.

# Load required modules
module load deeptools

# Name the region of interest file (bed) and signal files (bigWig). For
# each signal file, generate signal intensity matrix. The matrix will
# show signal intensity in a sample spaning 1000 base pairs before and
# after the center of each region of interest, in 10 base pair intervals.
regions="./outputs/SRR23310242.peaks.stringent.bed" 
bw_files=("./outputs/processed_data/SRR23310241.bw"
          "./outputs/processed_data/SRR23310242.bw")
for bw in ${bw_files[@]}; do
    out_name=$(echo ${bw} | sed "s/bw/matrix.tab/")
    computeMatrix reference-point -S ${bw} -R ${regions} -a 1000 \
    -b 1000 -bs 10 --referencePoint center --missingDataAsZero \
    --outFileName ${out_name}.gz --outFileNameMatrix ${out_name}
done
