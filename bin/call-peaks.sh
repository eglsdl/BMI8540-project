#!/bin/bash

#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10gb
#SBATCH --job-name=seacr
#SBATCH --error=./outputs/logs/seacr.err
#SBATCH --error=./outputs/logs/seacr_log.out

# This script is used to call peaks in the input intensity file
# to obtain bed file that lists regions of the genome that show
# high signal. In CUTAC-seq/ATAC-seq experiments, high signal
# corresponds to accessible chromatin. These accessible regions
# of interest are later used to plot the actual signal observed
# in different samples over them.

# Loading modules required by SEACR
module load bedtools/2.27
module load R/4.3

# Describing input file as well as path and prefix for output
# file.
in_file="./outputs/processed_data/SRR23310242.bg"
out_file="./outputs/SRR23310242.peaks"

# Calling peaks using SEACR. Arguments used:
# 0.01 - threshold of what fraction of peaks should be kept; IgG 
#        file could be used instead.
# non - when IgG file is available, it can be normalised to sample
#       before htrasholding. In this case using option none.
# stringent - take a more stringent approach in peak calling;
#             results in less peaks with smaller width.
./external/SEACR/SEACR_1.3.sh ${in_file} 0.01 non stringent ${out_file}
