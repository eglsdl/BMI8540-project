#!/bin/bash

#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10gb
#SBATCH --job-name=seacr
#SBATCH --chdir=/work/brahmalab/eglesidl/BMI8540-project
#SBATCH --error=/work/brahmalab/eglesidl/BMI8540-project/outputs/seacr.err
#SBATCH --output=/work/brahmalab/eglesidl/BMI8540-project/outputs/seacr_log.out

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
input_file="./outputs/mapped_data/SRR23310242.bg"
output_file="./outputs/SRR23310242.peaks"

# Calling peaks using SEACR. Arguments used:
# 0.01 - threshold of what fraction of peaks should be kept; IgG 
#        file could be used instead.
# non - when IgG file is available, it can be normalised to sample
#       before htrasholding. In this case using option none.
# stringent - take a more stringent approach in peak calling;
#             results in less peaks with smaller width.
./external/SEACR/SEACR_1.3.sh ${input_file} 0.01 non stringent ${output_file}
