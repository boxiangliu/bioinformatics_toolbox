#!/bin/bash 

# Usage: 
# 
# MODIFY: 
# -t 1-<num of jobs> 
# 
# CMD ARGS: 
# arg1	sorted.bam folder
# 
################# SCG settings ################### 
# Job Name 
#$ -N MarkDuplicate
# 
# Number of Jobs 
#$ -t 1-20
# tell the scheduler this takes a lot of memory 
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
# -cwd                                    
#
# set the output file
#$ -o  MarkDuplicate.log
#
# merge stdout and stderr                                         
#$ -j y
####################################################
 

set -u 
set -e 

# read command line arguments: 
input_dir=$1
output_dir=$1 
log=$output_dir/$(basename $0 .sh).log


# DO NOT MODIFY:
software=/srv/gs1/software/
hg19=/srv/gs1/projects/montgomery/bliu2/shared/ucsc_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
gatk=$software/gatk/gatk-3.3.0/GenomeAnalysisTK.jar
picard=$software/picard-tools/1.111 # e.g. $picard/MarkDuplicates.jar 
bwa=/srv/gs1/software/bwa/bwa-0.7.7/bin/bwa
star_genome=/srv/gs1/projects/montgomery/tnance/genomes/STAR/hg19_gencode14_overhang99
htseq=/srv/gs1/projects/montgomery/bliu2/tools/HTSeq-0.6.1/scripts
python=/home/bliu2/anaconda/bin/python
module load java 

# create array to store all sorted bam file names
cd $input_dir
inputs=(*.sorted.bam) # put the file extension here. 


# echo -e "$0 $*" | tee $log 
# echo "START: $(date)" | tee -a $log

i=$((SGE_TASK_ID-1))
output=${inputs[$i]/sorted.bam/dedup.bam} ; echo $output
metrics=${inputs[$i]/sorted.bam/metrics.txt}; echo $metrics 
java -Xmx4g -jar $picard/MarkDuplicates.jar INPUT=$output_dir/${inputs[$i]} OUTPUT=$output_dir/$output METRICS_FILE=$output_dir/$metrics
echo $output_dir/${inputs[$i]}.done 


# echo "FINISH: $(date)" | tee -a $log; touch template.done  

