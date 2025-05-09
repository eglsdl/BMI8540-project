#!/bin/bash

#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10gb
#SBATCH --job-name=plot
#SBATCH --error=./outputs/plot.err
#SBATCH --output=./outputs/plot_log.out

# This script calls python script plot-heatmaps.py to generate a plot
# for each sample of interest. Script takes matrix file and output
# file name as arguments.
mkdir ./outputs/plots
python3 ./bin/plot-heatmaps.py ./outputs/processed_data/SRR23310242.matrix.tab \
./outputs/plots/SRR23310242.png
python3 ./bin/plot-heatmaps.py ./outputs/processed_data/SRR23310241.matrix.tab \
./outputs/plots/SRR23310241.png
