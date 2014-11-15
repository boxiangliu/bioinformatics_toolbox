#!/bin/bash 

# Usage: 
# arg1	basecall directory
# arg2	alignment directory 

# set global variable:
log=$2/$(basename $0 .sh).log 
basecalls_dir=$1
alignments_dir=$2


# DO NOT MODIFY:

hg19=/srv/gs1/projects/montgomery/bliu2/shared/ucsc_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
BWAIndex=/srv/gs1/projects/montgomery/bliu2/shared/ucsc_hg19/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa
software=/srv/gs1/software/
gatk=$software/gatk/gatk-3.3.0/GenomeAnalysisTK.jar
picard=$software/picard-tools/1.111
bwa=/srv/gs1/software/bwa/bwa-0.7.7/bin/bwa


# main: 
cd $basecalls_dir
fastq_files=(*_L001_R1_001.fastq.trimmed.gz)

echo "START: $(date)"
echo "INPUT: $basecalls_dir" | tee $log 
echo "OUTPUT: $alignments_dir" | tee -a $log  

for fastq in ${fastq_files[*]}; do

	echo "START: $(date)" >> $log
	echo $fastq >> $log
	sam=${fastq/_L001_R1_001.fastq.trimmed.gz/.sam}
	$bwa mem -M -R '@RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit1' $BWAIndex $basecalls_dir/$fastq > $alignments_dir/$sam
	echo "FINISH: $(date)" >> $log 

done 