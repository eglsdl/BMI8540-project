#!/bin/bash

# This script is used to download files of interest from GSE224292.
# The input files are bigWig files generated from RNA Polymerase II 
# Serine 5 Phosphorylation CUTAC-seq experiment on mouse embryonic
# stem cells. 

# The data shows open chromatin that is associated with RNA Pol II
# S5P binding. The data is from two types of cells: 
# - 2i (cells in a more naïve pluripotent state);
# - SL (a mix of cells in pluripotent and differentiating state).

# Define directory to store files:
input_dir="./inputs"
cd ${input_dir}

# Download input files
echo "Downloading the sample file..."
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7019nnn/GSM7019043/suppl/GSM7019043%5FRNAPII%2DS5P%5FCUTAC%5F2i.mm10.bw
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7019nnn/GSM7019044/suppl/GSM7019044%5FRNAPII%2DS5P%5FCUTAC%5FSL.mm10.bw
echo "Done."
