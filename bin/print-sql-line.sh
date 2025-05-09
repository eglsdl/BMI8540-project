#!/bin/bash

# This script writes information about samples processed in this pipeline
# in a SQL format so that the information could be easily added to a
# database. If only one argument is provided and the argument is "new",
# the script will output code for creation of a new SQL database table.
# stored information is: SRR ID, total read count in raw fastq file,
# count of read pairs after trimming, count of successfully aligned reads,
# alignment rate (percentage), region of interest file name and final plot
# file name.
if [[ $# -eq 1 && $1 == "new" ]]; then
    echo "CREATE TABLE my_heatmap_records ("
    echo "    srr_id VARCHAR(11),"
    echo "    total_reads INT,"
    echo "    trimmed_reads INT,"
    echo "    aligned_reads INT,"
    echo "    alignment_rate VARCHAR(6),"
    echo "    peak_fn VARCHAR(100),"
    echo "    plot_fn VARCHAR(100),"
    echo "    date VARCHAR(10),"
    echo "    PRIMARY KEY(srr_id)"
    echo ");"
    exit 0
fi

# Otherwise, 3 arguments should be provided - SRR ID of sample processed,
# name of file with regions of interest used and name of the heatmap
# produced. The script will output an SQL code line to insert information
# about the sample into the database.
if [[ $# -eq 3 ]]; then
    srr=$1
    region_file=$2
    plot_file=$3
    total_reads=$(wc -l inputs/fasta/${srr}_1.fastq | awk '{print $1/4}')
    trimmed_reads=$(awk 'NR==1 {print $1}' \
	    ./outputs/processed_data/${srr}.bowtie2.out)
    aligned_reads=$(awk 'NR==4 || NR==5 {sum += $1} END {print sum}' \
	    ./outputs/processed_data/${srr}.bowtie2.out )
    aligned_prc=$(awk 'NR==6 {print $1}' \
	    ./outputs/processed_data/${srr}.bowtie2.out )
    date_added=$(date +"%Y-%m-%d")
    echo "INSERT INTO my_heatmap_records VALUES(\"${srr}\"," \
	    "${total_reads}, ${trimmed_reads}, ${aligned_reads}," \
	    "\"${aligned_prc}\", \"${region_file}\", \"${plot_file}\"," \
	    "\"${date_added}\");"
    exit 0
fi

# Otherwise, throw error that arguments provided are not correct.
if [[ ( $# -ne 1 || $1 != "new") && $# -ne 3 ]]; then
    echo "ERROR: Wrong arguments provided. To start new table use one" \
	    "argument \"new\". To add a new line to table, provide 3" \
	    "arguments: SRR ID, region BED file name and plot file name."
    exit 1
fi
