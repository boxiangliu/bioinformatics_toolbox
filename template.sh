#!/bin/bash 

# Usage: 
# arg1	
# arg2 

# read command line arguments: 
input_dir=$1
output_dir=$1 
log=$output_dir/$(basename $0 .sh)


# DO NOT MODIFY:
$software=/srv/gs1/software/
hg19=/srv/gs1/projects/montgomery/bliu2/shared/ucsc_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
gatk=$software/gatk/gatk-3.3.0/GenomeAnalysisTK.jar
picard=$software/picard-tools/1.111
bwa=/srv/gs1/software/bwa/bwa-0.7.7/bin/bwa


# create array to store all sorted bam file names
cd $input_dir
inputs=(*.) # put the file extension here. 


echo -e $* | tee $log 
echo "START: $(date)" | tee -a $log

for input in ${inputs[*]}; do 
	# your code here: 

done 

echo "FINISH: $(date)" | tee -a $log 

