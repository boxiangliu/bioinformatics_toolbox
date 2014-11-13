#!/bin/bash 
set -u 
set -e 

# first command line input is the input directory
input_dir=$1 

# second is output directory
output_dir=$2

# create an array of all fastq files:  
fastq_files=(*_L001_R1_001.fastq.gz)
length=${#fastq_files[*]}
echo "$length files..."
#### call fastqc #### 
module load fastqc/0.11.2

for fastq in $(fastq_files); do 
	echo $fastq 
	fastqc $input_dir/$fastq --outdir=$output_dir
done 


