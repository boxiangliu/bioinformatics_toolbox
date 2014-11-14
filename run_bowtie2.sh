#!/bin/bash
# INPUT
# arg1	directory where fastq files are located 
# arg2	directory where sam files will be written


#### scg3 parameters #### 
# -t 1
#
# set the maximum memory usage (per slot)
# -l h_vmem=8G
#
# set the number of cores
# -pe shm 8
#
# set the maximum run time
#$ -l h_rt=06:00:00
# 
# check for errors in the job submission options
#$ -w e
# merge the standard error stream into the starndard output stream 
# -j y
# set the path for standard output stream 
#$ -o  /srv/gs1/projects/montgomery/bliu2/ancestry/log
# check for errors in the job submission options
#$ -w e
#########################

# load bowtie2: 
module load bowtie/2.2.1

# create directory aliases: 
fastq_dir=$1
sam_dir=$2
index_dir=/srv/gs1/projects/montgomery/bliu2/shared/bowtie2 # change this to your bowtie2 index directory.
log_dir=$2
echo "INPUT: $fastq_dir" > $sam_dir/run_bowtie2.log
echo "OUTPUT: $sam_dir" >> $sam_dir/run_bowtie2.log

# create an array of fastq files: 
cd $fastq_dir
fastq_files=(*_L001_R1_001.fastq.gz)
length=${#fastq_files[*]}

echo "$length files..."
echo "START: $(date)"  

# for file in 
# bowtie2 -x $index_dir/hg19 -U $fastq_dir/${fastq[$SGE_TASK_ID]} -S $sam_dir/${fastq[$SGE_TASK_ID]}.sam > $log_dir/${fastq[$SGE_TASK_ID]}.log
# echo 'Done.'

for fastq in ${fastq_files[*]}; do 
	
	sam=${fastq/_L001_R1_001.fastq.gz/.sam}
	log=${fastq/_L001_R1_001.fastq.gz/.log}

	echo "START: $(date)" >> $sam_dir/run_bowtie2.log
	echo $fastq >> $sam_dir/run_bowtie2.log

	bowtie2 -x $index_dir/hg19 -U $fastq_dir/$fastq -S $sam_dir/$sam 2>$log_dir/$log

	echo "FINISH: $(date)" >> $sam_dir/run_bowtie2.log

done

echo "FINISH: $(date)"


