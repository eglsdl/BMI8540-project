#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=20gb
#SBATCH --job-name=download-process
#SBATCH --chdir=/work/brahmalab/eglesidl/BMI8540-project
#SBATCH --error=/work/brahmalab/eglesidl/BMI8540-project/outputs/pre-process.err
#SBATCH --output=/work/brahmalab/eglesidl/BMI8540-project/outputs/pre-process_log.out

# This script is used to download and pre-process the files of interest on UNL HCC Swan.

# These input files are sequencing data from CUTAC-seq experiment targeting
# RNA Polymerase II Serine 5 Phosphorylation on mouse embryonic stem cells.
# The cells are from SL cells, meaning that a mix of cells in pluripotent
# and differentiating state can be found in each sample.
# SRR23310241 - 8 hour Flavopiridol (RNA Pol II inhibitor) treatment.
# SRR23310242 - no treatment.

# Load required tools
module load SRAtoolkit/2.11
module load cutadapt
module load bowtie/2.5
module load biodata

# Define directory to store input and output files:
fasta_input_dir="./inputs/fasta"
cutadapt_output_dir="./outputs/cutadapt"
bowtie2_output_dir="./outputs/bowtie2"
mkdir ${fasta_input_dir} ${cutadapt_output_dir} ${bowtie2_output_dir}

# Add SRR IDs of my files of interest into a list and loop thorugh
srr_ids=("SRR23310241" "SRR23310242")
for srr in ${srr_ids[@]}; do

    # Download input files using SRA-toolkit
    echo "Downloading FASTA files for ${srr}..."
    prefetch ${srr} -O ${fasta_input_dir}
    fasterq-dump --split-files ${srr} -O ${fasta_input_dir} --threads 8
    echo "Done."

    # Trim adapters using cutadapt
    echo "Trimming adapters for ${srr}..."
    fq1="${fasta_input_dir}/${srr}_1.fastq"
    fq2="${fasta_input_dir}/${srr}_2.fastq"
    trim1="${cutadapt_output_dir}/${srr}_1.trimmed.fastq"
    trim2="${cutadapt_output_dir}/${srr}_2.trimmed.fastq"
    cutadapt -j 8 -m 20 --nextseq-trim 20 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -Z -o ${trim1} -p ${trim2} ${fq1} ${fq2} \
    >& ${cutadapt_output_dir}/${srr}.cutadapt
    echo "Done."

    # Map reads to Mus Musculus (mm10) genome using Bowtie2
    echo "Mapping ${srr} to Mus Musculus..."
    sam_mm10="${bowtie2_output_dir}/${srr}.mm10.sam"
    out_mm10="${bowtie2_output_dir}/${srr}.mm10.bowtie2.out"
    bowtie2 --very-sensitive-local --soft-clipped-unmapped-tlen \
    --no-unal --no-mixed --no-discordant --dovetail --phred33 -I 10 -X 1000 \
    --threads 8 -x $BOWTIE2_MUS_MUSCULUS_UCSC_MM10 -1 ${trim1} -2 ${trim2} \
    > ${sam_mm10} 2> ${out_mm10}
    echo "Done."

    # Map reads to Drosophila Melanogaster genome using Bowtie2
    echo "Mapping ${srr} to Drosophila Melanogaster..."
    sam_dros="${bowtie2_output_dir}/${srr}.dros.sam"
    out_dros="${bowtie2_output_dir}/${srr}.dros.bowtie2.out"
    bowtie2 --very-sensitive-local --soft-clipped-unmapped-tlen \
    --no-unal --no-mixed --no-discordant --dovetail --phred33 -I 10 -X 1000 \
    --threads 8 -x $BOWTIE2_DROSOPHILA_MELANOGASTER_NCBI_BUILD5_41 -1 \
    ${trim1} -2 ${trim2} > ${sam_dros} 2> ${out_dros}
    echo "Done."
done
