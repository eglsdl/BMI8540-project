#!/bin/bash

#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10gb
#SBATCH --job-name=dwn-prc
#SBATCH --error=./outputs/logs/pre-process.err
#SBATCH --output=./outputs/logs/pre-process_log.out

# This script is used to download data of interest, trim adapters,
# map it to Mus Musculus (mm10) genome, convert mapping results to
# bed file format, read depth normalise samples and obtain bedGraph
# and bigWig files that show signal strength. Settings for adapter 
# trimming and mapping were mirrored from instructions from GEO:
# SRR23310241 - GSM7019046
# SRR23310242 - GSM7019045 

# Noting the job start time in log
echo "Job started at:"
date +"%Y-%m-%d %H:%M:%S"

# Loading required tools
module load bedtools
module load biodata
module load bowtie/2.4
module load cutadapt
module load samtools
module load SRAtoolkit/2.11
module load ucsc-bedgraphtobigwig/455

# Defining directories to store input and output files and create them:
fasta_input_dir="./inputs/fasta"
output_dir="./outputs/processed_data"
mkdir -p ${fasta_input_dir} ${output_dir} ${output_dir}

# For mapping, pre-built Mus Musculus UCSC genome (mm10) will be used. It
# is stored on HCC Swan biodata module and can be accessed using variable
# $BOWTIE2_MUS_MUSCULUS_UCSC_MM10. Chromosome and genome length files
# can be obtained with code lines bellow. To  generate the chr_length.tab
# file I am using awk to loop through lines and count base pairs (bp) in
# each chromosome:
# -- If a line starts with symbol ">" it's a new chomosome - name of
# previous chromosome and it's bp count is printed; name of the new
# chromosome is saved;
# -- If a line does not start with ">", it's length is added to a counter;
# -- Afterwards, the information is sorted based on chromosome name. 
chr_length="./inputs/chr_length.tab"
bowtie2-inspect $BOWTIE2_MUS_MUSCULUS_UCSC_MM10 | awk '/^>/ { \
	if (chr != "") print chr "\t" seqlen; chr=substr($0,2); \
	seqlen=0; next} {seqlen += length($0)} END \
       	{print chr "\t" seqlen}' | sort -k1,1V > ${chr_length}
genome_length=$(awk '{sum += $2} END {print sum}' ${chr_length})

# Adding SRR IDs of the files of interest into a list and loop through.
# These input files are sequencing data from a CUTAC experiment targeting
# RNA Polymerase II Ser 5 Phosphorylation on mouse embryonic stem cells.
# SRR23310241 - 8 hour Flavopiridol (RNA Pol II inhibitor) treatment.
# SRR23310242 - no treatment.
srr_ids=("SRR23310241" "SRR23310242")
for srr in ${srr_ids[@]}; do

    # Downloading input files using SRA-toolkit
    echo "Downloading FASTA files for ${srr}..."
    prefetch ${srr} -O ${fasta_input_dir}
    fasterq-dump --split-files ${srr} -O ${fasta_input_dir} --threads 4
    echo "Done."

    # Trimming adapters using cutadapt. Settings were mirrored from
    # instructions from GEO (GSM7019046 & GSM7019045)
    echo "Trimming adapters for ${srr}..."
    fq1="${fasta_input_dir}/${srr}_1.fastq"
    fq2="${fasta_input_dir}/${srr}_2.fastq"
    trim1="${output_dir}/${srr}_1.trimmed.fastq"
    trim2="${output_dir}/${srr}_2.trimmed.fastq"
    cutadapt -j 4 -m 20 --nextseq-trim 20 \
	    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
	    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
	    -Z -o ${trim1} -p ${trim2} ${fq1} ${fq2} \
	    >& ${output_dir}/${srr}.cutadapt
    echo "Done."

    # Mapping reads to Mus Musculus (mm10) genome using Bowtie2. Settings
    # were mirrored from instructions from GEO (GSM7019046 & GSM7019045)
    echo "Mapping ${srr} to mm10..."
    sam="${output_dir}/${srr}.sam"
    out="${output_dir}/${srr}.bowtie2.out"
    bowtie2 --very-sensitive-local --soft-clipped-unmapped-tlen \
	    --no-mixed --no-discordant --dovetail -q --phred33 -I 10 \
	    -X 1000 --threads 4 -x $BOWTIE2_MUS_MUSCULUS_UCSC_MM10 \
	    -1 ${trim1} -2 ${trim2} > ${sam} 2> ${out}
    echo "Done."

    # Extracting the properly paired alignments into a bed file. To do
    # that file is first converted into bedpe format (both paired end (pe)
    # reads are output in the same line). Then, using cut, chromosome
    # name, start point of one pe read, end point of other pe read and
    # gene information is output. The file is sorted and fragment length
    # information is added to the 5th column using awk.
    echo "Converting mapped ${srr} sam file to bed..."
    bed="${output_dir}/${srr}.bed"
    samtools view -@ 4 -b -h ${sam} | bedtools bamtobed -bedpe -i stdin |\
	    cut -f1,2,6,7 | sort -k1,1 -k2,2n -k3,3n | awk -v OFS='\t' \
	    '{len = $3 - $2 ; print $0, len }' > ${bed}
    echo "Done."

    # Keeping only reads that are 10-1000 bp long. Calculating a 
    # scale factor for read depth normalisation by dividing total
    # bp count in mm10 genome by count of mapped reads.
    echo "Generating scale factor for ${srr}..."
    bed_lim="${output_dir}/${srr}.10-1000.bed"
    awk '{if ($5 >= 10 && $5 <= 1000) print }' ${bed} > ${bed_lim}
    count=$(awk '{sum += $5} END {print sum}' ${bed_lim})
    scale_factor=$(echo "${genome_length}/${count}" | bc -l)
    echo "Scale factor: ${scale_factor}."
    echo "Mapped read count: ${count}."
    echo "Done."

    # Extracting signal intensity information from the samples and
    # normalising it to read depth. The bedGraph file is then converted
    # to a bigWig file.
    echo "Generating bedGraph and bigWig files for ${srr}..."
    bg="${output_dir}/${srr}.bg"
    bw="${output_dir}/${srr}.bw"
    bedtools genomecov -bg -scale ${scale_factor} -i ${bed_lim} \
	    -g ${chr_length} > ${bg}
    bedGraphToBigWig ${bg} ${chr_length} ${bw}
    echo "Done." 
done

# Noting job finish time in logs
echo "Job finished at:" 
date +"%Y-%m-%d %H:%M:%S"
