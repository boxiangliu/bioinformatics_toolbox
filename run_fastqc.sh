#!/bin/bash 

# Usage:
# arg1	fastq directory 
# arg2	fastQC directory
# arg3	filename, e.g. "*.fastq.gz". MAKE SURE TO INCLUE THE QUOTATION MARK "". 


set -u 
set -e 


# first command line input is the input directory
input_dir=$1 

# second is output directory
output_dir=$2
log=$output_dir/$(basename $0 .sh).log

# create an array of all fastq files: 
cd $input_dir 
fastq_files=($3)
length=${#fastq_files[*]}
echo "$length files..."
#### call fastqc #### 
module load fastqc/0.11.2

echo -e $* | tee $log 
echo "START: $(date)" | tee -a $log

for fastq in ${fastq_files[*]}; do 
	fastqc $input_dir/$fastq --outdir=$output_dir | tee -a $log
done 


echo "FINISH: $(date)" | tee -a $log 
