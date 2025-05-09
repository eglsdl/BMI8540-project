#!/bin/bash

#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10gb
#SBATCH --job-name=sql
#SBATCH --error=./outputs/sql.err
#SBATCH --output=./outputs/sql_log.out

# This script calls bash script print-sql-line.sh to generate a SQL
# script for each processed sample. The SQL script can later be used
# to create a database table and store relevant data. Script takes
# keyword "new" as input to generate new database table or 3 arguments:
# SRR ID, region of interest file name and output plot file name, to
# produce insert statements for the sample.
output_file="./outputs/store_sample_info.sql"
./bin/print-sql-line.sh new > ${output_file}
./bin/print-sql-line.sh SRR23310241 outputs/SRR23310242.peaks.stringent.bed outputs/plots/SRR23310241.png >> ${output_file}
./bin/print-sql-line.sh SRR23310242 outputs/SRR23310242.peaks.stringent.bed outputs/plots/SRR23310242.png >> ${output_file}
echo "File ${output_file} generated."
