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
input_extension="sam"
output_extension="sorted.bam"

################# SCG settings ################### 
# Job Name 
#$ -N compress_sort_and_index 
# 
# Array Job 
#$ -t 1-97
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


# input_dir=$1
# output_dir=$2 
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
inputs=(*$input_extension) # put the file extension here. 
i=$((SGE_TASK_ID-1))
sam={inputs[$i]/$input_extension/$output_extension}
output=${inputs[$i]/$input_extension/$output_extension}

start=$(date)

java -Xmx2g -jar $picard/SortSam.jar INPUT=$input_dir/${inputs[$i]} OUTPUT=/dev/stdout SORT_ORDER=coordinate | \
java -Xmx2g -jar $picard/SamFormatConverter.jar INPUT=/dev/stdin OUTPUT=$output_dir/$output  CREATE_INDEX=true

finish=$(date)
touch $output_dir/$output.done; echo -e "START: $start\nFINISH: $finish" > $output_dir/$output.done