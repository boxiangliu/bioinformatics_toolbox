#!/bin/bash 
# INPUT: 
# arg1	input directory with fastq files 
# arg2	output directory where you want the trimmed fastq files to go. 


set -u 
set -e 

# first command line input is the input directory
input_dir=$1 

# second is output directory
output_dir=$2

# adapter sequences, default to forward-read TruSeq adapter: 
adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG"

# create an array of all fastq files:  
cd $input_dir
fastq_files=(*_L001_R1_001.fastq.gz)
length=${#fastq_files[*]}

echo "File List\n" >> $output_dir/run_cutadapt.log 
echo ${fastq_files[*]} >> $output_dir/run_cutadapt.log 
echo "$length files..."

# invoke cutadapt over all files in array: 
for orig in ${fastq_files[*]}; do 
	trimmed=${orig/fastq.gz/fastq.trimmed.gz}
	log=${orig/fastq.gz/log}
	cutadapt -a $adapter -o $output_dir/$trimmed $input_dir/$orig > $output_dir/$log
done 


