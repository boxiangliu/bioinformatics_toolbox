#!/bin/sh 
#
# Usage: 
# 
# SCG SETTINGS:
# -t 1:<num of jobs>
# 
# CMD ARGS
# arg1	fastq directory
# arg2	trimmed fastq directory 
# arg3	fastqc directory, default to arg2
# arg4	s (single) or p (paired-end), default to p. 
# 
# MODIFY:
read1_extension=_R1_001.fastq.gz
read2_extension=_R2_001.fastq.gz
min_qual=30

################# SCG settings ################### 
# Job Name 
#$ -N run_trim_galore
# 
# Array Job 
#$ -t 1-51
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
# working directory      
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

module load trim_galore/0.3.3

# read CMD ARGS: 
input_dir=$1
output_dir=$2
fastqc_dir=${3:-$output_dir}
read_mode=${4:-"p"}
log=$output_dir/$(basename $0 .sh).log

# log commnad and time: 
echo -e "COMMAND: $0 $*" | tee $log
echo "START: $(date)" | tee -a $log

# modify the suffix for R1 and R2:
cd $input_dir
fastq_R1=(*$read1_extension)
fastq_R2=(*$read2_extension)

# adapter sequences, default to forward-read TruSeq adapter: 
adapter_full="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG"
adapter="AGATCGGAAGAGC"
adapter2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
adapter2_mmPCR="TTACTATGCCGCTGGTGGCTGTGAGAAAGGGATGTGCTCGCAATAGCTCCAG" # mmPCR uses a customized forward adapter. 

# main: 
i=$((SGE_TASK_ID-1))
read1=${fastq_R1[$i]}
read2=${read1/R1/R2}
prefix=$(basename $read1 $read1_extension)

if [ $read_mode == "p" ]; then 
	echo "Paired-ended reads." 
	trim_galore --paired -q $min_qual --phred33 --fastqc --fastqc_args "--output_dir $fastqc_dir" --output_dir $output_dir --adapter $adapter --adapter2 $adapter2 --stringency 1 $read1 $read2
	# --paired 		paired-end reads
	# -q			quality cutoff for trimming, default 20
	# --phred33		Sanger/Illumina 1.9+ encoding 
	# --fastqc_dir 	run fastqc after trimming
	# --fastqc_args fastqc arguments 
	# --adapter 	adapter sequences 
	# --stringency 	minimum number of adapter-overlapping bases for trimming

	
elif [ $read_mode == "s" ]; then
	echo "Single-ended reads." 
	trim_galore -q $min_qual --phred33 --fastqc --fastqc_args "--output_dir $fastqc_dir" --output_dir $output_dir --adapter $adapter --stringency 1 $read1
	# -q			quality cutoff for trimming, default 20
	# --phred33		Sanger/Illumina 1.9+ encoding 
	# --fastqc_dir 	run fastqc after trimming
	# --fastqc_args fastqc arguments 
	# --adapter 	adapter sequences 
	# --stringency 	minimum number of adapter-overlapping bases for trimming
	
else 
	echo "Run mode is neither [s]ingle or [p]aired-end." ; exit 1

fi

echo "FINISH: $(date)" | tee -a $log ; touch $prefix.done