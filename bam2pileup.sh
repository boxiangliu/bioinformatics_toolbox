#!/bin/bash

# Useage: 
# This script convert all bam files in a folder to pileup format. 
# 
# INPUTS 
# arg1	bams directory. It's okay if there are files other than bam. The script will look for *.sorted.bam extension. 
# arg2 	pileup directory  
# arg3	target sites file 

module load samtools/0.1.19

alignments_dir=$1
pileup_dir=$2
target_sites_filename=${3:-"None"}
hg19=/srv/gs1/projects/montgomery/shared/genome/hg19/hg19.fa
log=$pileup_dir/bam2pileup.log 

echo "START: $(date)"

echo "INPUT: $alignments_dir" | tee $log 
echo "OUTPUT: $pileup_dir" | tee -a $log 

cd $alignments_dir
bam_files=(*.sorted.bam)
echo "${#bam_files[*]} files..." | tee -a $log  # number of files 

if [ $target_sites_filename != "None" ]; then
	echo "Pileup over target sites in $target_sites_filename..."
	for bam_file in ${bam_files[*]}; do 

		echo "START: $(date)" >> $log
		echo $bam_file >> $log
		pileup=${bam_file/sorted.bam/pileup} # get filename without extension 
		samtools mpileup -B -f $hg19 -l $target_sites_filename $alignments_dir/$bam_file > $pileup_dir/$pileup
		echo "FINISH: $(date)" >> $log 

	done
else 
	echo "No target sites file specified. Pileup over all sites..."
	for bam_file in ${bam_files[*]}; do 

		echo "START: $(date)" >> $log
		echo $bam_file >> $log
		pileup=${bam_file/sorted.bam/pileup} # get filename without extension 
		samtools mpileup -B -f $hg19 $alignments_dir/$bam_file > $pileup_dir/$pileup
		echo "FINISH: $(date)" >> $log 

	done
fi 

echo "FINISH: $(date)" | tee -a $log; touch $pileup.done 