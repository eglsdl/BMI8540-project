#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10gb
#SBATCH --job-name=download-process
#SBATCH --chdir=/work/brahmalab/eglesidl/BMI8540-project
#SBATCH --error=/work/brahmalab/eglesidl/BMI8540-project/outputs/pre-process.err
#SBATCH --output=/work/brahmalab/eglesidl/BMI8540-project/outputs/pre-process_log.out

# This script is used to download and pre-process the files of interest.

# These input files are sequencing data from CUTAC-seq experiment targeting
# RNA Polymerase II Serine 5 Phosphorylation on mouse embryonic stem cells.
# The cells are from SL cells, meaning that a mix of cells in pluripotent
# and differentiating state can be found in each sample.
# SRR23310241 - 8 hour Flavopiridol (RNA Pol II inhibitor) treatment.
# SRR23310243 - no treatment.

# Loading SRA-toolkit to use for data download
module load SRAtoolkit/2.11

# Define directory to store files:
fasta_input_dir="./inputs/fasta"
mkdir ${fasta_input_dir}

# Download input files
echo "Downloading sample files..."
prefetch SRR23310241 -O ${fasta_input_dir}
fasterq-dump --split-files SRR23310241 -O ${fasta_input_dir} --threads 4
prefetch SRR23310243 -O ${fasta_input_dir}
fasterq-dump --split-files SRR23310243 -O ${fasta_input_dir} --threads 4
#wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7019nnn/GSM7019043/suppl/GSM7019043%5FRNAPII%2DS5P%5FCUTAC%5F2i.mm10.bw
#wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7019nnn/GSM7019044/suppl/GSM7019044%5FRNAPII%2DS5P%5FCUTAC%5FSL.mm10.bw
echo "Done."
