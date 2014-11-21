#!/bin/bash 

# Usage: 
# arg1	
# arg2 

set -u 
set -e 

# user parameters:
type=sam
readlength=76

# read command line arguments: 
input_dir=$1
output_dir=$1 
log=$output_dir/$(basename $0 .sh).log


# DO NOT MODIFY:
software=/srv/gs1/software/
hg19=/srv/gs1/projects/montgomery/bliu2/shared/ucsc_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
gatk=$software/gatk/gatk-3.3.0/GenomeAnalysisTK.jar
picard=$software/picard-tools/1.111
bwa=/srv/gs1/software/bwa/bwa-0.7.7/bin/bwa
star_genome=/srv/gs1/projects/montgomery/tnance/genomes/STAR/hg19_gencode14_overhang99
htseq=/srv/gs1/projects/montgomery/bliu2/tools/HTSeq-0.6.1/scripts
python=/home/bliu2/anaconda/bin/python

# create array to store all sorted bam file names
cd $input_dir
inputs=(*.sam) # put the file extension here. 


echo -e "$0 $*" | tee $log 
echo "START: $(date)" | tee -a $log

for input in ${inputs[*]}; do 

	$python $htseq/htseq-qa -t $type -r $readlength $input | tee -a $log 

done 

echo "FINISH: $(date)" | tee -a $log; touch run_htseq_qa.done 

