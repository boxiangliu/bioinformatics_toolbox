#!/bin/bash 

# Usage: 
# 
# Run qsub array jobs on all files in a directory. 
# To uses this script, 
# 1. change the second parameters of -t to the number of jobs you need to run. 
# 2. assign your input file extension to the variable extension, and output file extension to output_extension.
# 3. add you command on line "your comment here:"
# 4. if your array jobs require large memory, enable -P large_mem and -l h_vmem=40G.
# 
# SCG SETTINGS: 
# -t 1:<num of jobx>
# 
# CMD ARGS:
# -i	input dir
# -o	output dir
# 
# MODIFY: 
extension="_trimmed.fq.gz"
output_extension=".sam"

################# SCG settings ################### 
# Job Name 
#$ -N BWA_MEM
# 
# Array Job 
#$ -t 1-20
# 
# Request Large Memory Machine  
# -P large_mem
# 
# memory usage 
# -l h_vmem=40G 
#
# maximum run time 
#$ -l h_rt=96:00:00  
# 
# check for errors in the job submission options
#$ -w e
# 
# run on multiple threads                     
#$ -pe shm 4                                 
#
# run job in current working directory      
#$ -cwd                                    
#
# set the output file
# -o  template.log
#
# merge stdout and stderr                                         
#$ -j y
####################################################


set -u 
set -e 

# read command line arguments: 
parrot=true
while getopts ":i:o:" opt; do 
	case $opt in 
		i) input_dir=$OPTARG
		[ $parrot ] & echo "INPUT: $input_dir"
		;;
		o) output_dir=$OPTARG
		[ $parrot ] & echo "OUTPUT: $output_dir"
		;;
		\?) echo "Invalid option -$OPTARG"; exit 1
		;;
	esac
done 

log=$output_dir/$(basename $0 .sh).log


# DO NOT MODIFY:
software=/srv/gs1/software/
hg19=/srv/gs1/projects/montgomery/bliu2/shared/ucsc_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
gatk=$software/gatk/gatk-3.3.0/GenomeAnalysisTK.jar
picard=$software/picard-tools/1.111
bwa=/srv/gs1/software/bwa/bwa-0.7.7/bin/bwa
BWAIndex=/srv/gs1/projects/montgomery/bliu2/shared/ucsc_hg19/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa
star_genome=/srv/gs1/projects/montgomery/tnance/genomes/STAR/hg19_gencode14_overhang99
htseq=/srv/gs1/projects/montgomery/bliu2/tools/HTSeq-0.6.1/scripts
python=/home/bliu2/anaconda/bin/python
gencode14=/srv/gs1/projects/montgomery/shared/annotation/gencode.v14.annotation.gtf
gencode21=/srv/gs1/projects/montgomery/shared/annotation/gencode.v21.annotation.gtf
module load java 
bt=/srv/gs1/projects/montgomery/bliu2/bioinformatics_toolbox
brt=/srv/gs1/projects/montgomery/bliu2/brt

# create array to store all sorted bam file names
cd $input_dir
inputs=(*$extension) # put the file extension here. 
i=$((SGE_TASK_ID-1))
output=${inputs[$i]/$extension/$output_extension}

start=$(date)
#
# Using hg19:
#
# $bwa mem -t 8 -M -R '@RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit1' $BWAIndex $input_dir/${inputs[$i]} > $output_dir/$output
#
# Using custom genome: 
#
$bwa mem -t 8 -M -R '@RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit1' /srv/gs1/projects/montgomery/bliu2/ancestry/data/Bosh_ASW_mmPCR_2/targeted_reference_genome/hg19_150_aims.fa $input_dir/${inputs[$i]} > $output_dir/$output

finish=$(date)

touch $output_dir/$output.done; echo -e "START: $start\nFINISH: $finish" > $output_dir/$output.done



############## qlogin mode, for debugging #################



#!/bin/bash 

# Usage: 
# arg1	basecall directory
# arg2	alignment directory 

# set global variable:
# log=$2/$(basename $0 .sh).log 
# basecalls_dir=$1
# alignments_dir=$2


# # DO NOT MODIFY:

# hg19=/srv/gs1/projects/montgomery/bliu2/shared/ucsc_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
# BWAIndex=/srv/gs1/projects/montgomery/bliu2/shared/ucsc_hg19/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa
# software=/srv/gs1/software/
# gatk=$software/gatk/gatk-3.3.0/GenomeAnalysisTK.jar
# picard=$software/picard-tools/1.111
# bwa=/srv/gs1/software/bwa/bwa-0.7.7/bin/bwa


# # main: 
# cd $basecalls_dir
# fastq_files=(*_L001_R1_001.fastq.trimmed.gz)

# echo "START: $(date)"
# echo "INPUT: $basecalls_dir" | tee $log 
# echo "OUTPUT: $alignments_dir" | tee -a $log  

# for fastq in ${fastq_files[*]}; do

# 	echo "START: $(date)" >> $log
# 	echo $fastq >> $log
# 	sam=${fastq/_L001_R1_001.fastq.trimmed.gz/.sam}
# 	$bwa mem -M -R '@RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit1' $BWAIndex $basecalls_dir/$fastq > $alignments_dir/$sam
# 	echo "FINISH: $(date)" >> $log 

# done 