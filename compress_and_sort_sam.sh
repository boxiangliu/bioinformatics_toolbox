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
extension="sam"
output_extension="sorted"

################# SCG settings ################### 
# Job Name 
#$ -N compress_and_sort
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
star_genome=/srv/gs1/projects/montgomery/tnance/genomes/STAR/hg19_gencode14_overhang99
htseq=/srv/gs1/projects/montgomery/bliu2/tools/HTSeq-0.6.1/scripts
python=/home/bliu2/anaconda/bin/python
gencode14=/srv/gs1/projects/montgomery/shared/annotation/gencode.v14.annotation.gtf
gencode21=/srv/gs1/projects/montgomery/shared/annotation/gencode.v21.annotation.gtf
module load java 
bt=/srv/gs1/projects/montgomery/bliu2/bioinformatics_toolbox
brt=/srv/gs1/projects/montgomery/bliu2/brt
module load samtools/0.1.19

# create array to store all sorted bam file names
cd $input_dir
inputs=(*$extension) # put the file extension here. 
i=$((SGE_TASK_ID-1))
output=${inputs[$i]/$extension/$output_extension}

# your command here: 
samtools view -bS $input_dir/${inputs[$i]} | samtools sort - $output_dir/$output


touch $output_dir/$output.done  



################## Below is the qlogin version for debugging ####################

# INPUT
# arg1	sam directory. It is okay if other files are also in the folder. The script will look for *.sam exntension. 


# sam_dir=$1
# bam_dir=$1 
# echo "INPUT: $sam_dir" > $sam_dir/compress_and_sort_sam.log
# echo "START: $(date)"

# cd $sam_dir
# sam_files=(*.sam)
# echo ${#sam_files[*]}
# echo ${sam_files[*]}

# for sam_file in ${sam_files[*]}; do 
# 	echo "START: $(date)" >> $sam_dir/compress_and_sort_sam.log
# 	echo "INPUT: $sam_file" >> $sam_dir/compress_and_sort_sam.log

#  	bam_file=${sam_file/.sam/.sorted} # get filename without extension 

#  	echo "OUTPUT: $bam_file" >> $sam_dir/compress_and_sort_sam.log

# 	samtools view -bS $sam_dir/$sam_file | samtools sort - $sam_dir/$bam_file
	
# 	echo -e "FINISH: $(date)\n" >> $sam_dir/compress_and_sort_sam.log

# done

# echo "FINISH: $(date)"
# touch $bam_dir/compress_and_sort_sam.done