#!/bin/bash

# INPUTS 
# arg1	a directory where sorted bams files are located. It's okay if there are files other than bam. The script will look for *.sorted.bam extension. 
# arg2 	a directory where you want to put output bcf files. 
module load samtools/0.1.19

alignments_dir=$1
variants_dir=$2
hg19=/srv/gs1/projects/montgomery/shared/genome/hg19/hg19.fa
log=$variants_dir/call_genotype_likelihood.log 

echo "START: $(date)"

echo "INPUT: $alignments_dir" | tee $log 
echo "OUTPUT: $variants_dir" | tee -a $log 

cd $alignments_dir
bam_files=(*.sorted.bam)
echo "${#bam_files[*]} files..." | tee -a $log  # number of files 

for bam_file in ${bam_files[*]}; do 

	echo "START: $(date)" >> $log
	echo $bam_file >> $log
	bcf=${bam_file/sorted.bam/bcf} # get filename without extension 
	samtools mpileup -g -f $hg19 $alignments_dir/$bam_file > $variants_dir/$bcf
	echo "FINISH: $(date)" >> $log 

done

echo "FINISH: $(date)" | tee -a $log 