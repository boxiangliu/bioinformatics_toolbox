#!/bin/sh 
#
# Usage: 
# 
# SCG SETTINGS:
# -t 1:<num of jobs>
# 
# CMD ARGS
# arg1	fastq directory
# arg2	alignment directory 
# arg3	s (single) or p (paired-end), default to p. 
# 
# MODIFY:
read1_extension=_R1_001.fastq.gz
read2_extension=_R2_001.fastq.gz


################# SCG settings ################### 
# Job Name 
#$ -N STAR
# 
# Array Job 
#$ -t 1-51
# 
# Request Large Memory Machine  
#$ -P large_mem
# 
# memory usage 
#$ -l h_vmem=40G 
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

# read command line arguments: 
input_dir=$1
cd $input_dir
output_dir=$2
log=$output_dir/$(basename $0 .sh).log
read_mode=${3:-"p"}

# modify the suffix for R1 and R2:
fastq_R1=(*$read1_extension)
fastq_R2=(*$read2_extension)

# DO NOT MODIFY:
software=/srv/gs1/software/
hg19=/srv/gs1/projects/montgomery/bliu2/shared/ucsc_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
gatk=$software/gatk/gatk-3.3.0/GenomeAnalysisTK.jar
picard=$software/picard-tools/1.111
bwa=/srv/gs1/software/bwa/bwa-0.7.7/bin/bwa
star_genome_v14=/srv/gs1/projects/montgomery/tnance/genomes/STAR/hg19_gencode14_overhang99
star_genome_v21=/srv/gs1/projects/montgomery/bliu2/shared/ucsc_hg19/Homo_sapiens/UCSC/hg19/Sequence/STARGenome/hg19_gencode21_overhang75
module load STAR/2.3.0e

#### main: 

# check the number of forward and reverse files match:
num_R1=${#fastq_R1[*]}
num_R2=${#fastq_R2[*]}
# [ $num_R1 -eq $num_R2 ] || (echo "Number of forward read files and reverse read files not equal" > &2 ; exit 1)

echo -e "COMMAND: $0 $*" | tee $log
echo "START: $(date)" | tee -a $log

# load genome into shared memory:
STAR --genomeDir $star_genome_v21 --genomeLoad LoadAndExit; echo "FINISH LOADING GENOME $(date)" | tee -a $log 

## RUNNING SEQUENTIAL JOBS;
# for read1 in ${fastq_R1[*]}; do 
# 	read2=${read1/R1/R2}
# 	prefix=$(basename $read1 _R1_001.fastq.gz)
# 	echo -e "SAMPLES: $read1\t$read2" 
# 	echo "OUTPUT: $output_dir/${prefix}"
# 	STAR --genomeDir $star_genome --readFilesIn $input_dir/$read1 $input_dir/$read2 --readFilesCommand zcat --outFileNamePrefix $output_dir/${prefix}. --runThreadN 10
# done

# RUNNING PARALLEL/ARRAY JOBS: 
i=$((SGE_TASK_ID-1))
read1=${fastq_R1[$i]}
read2=${read1/R1/R2}
prefix=$(basename $read1 $read1_extension)

if [ $read_mode == "q" ]; then 
	STAR --genomeDir $star_genome --readFilesIn $input_dir/$read1 $input_dir/$read2 --readFilesCommand zcat --outFileNamePrefix $output_dir/${prefix}. --runThreadN 10
	echo "Aligning single-ended reads." 
elif [ $read_mode == "s" ]; then
	STAR --genomeDir $star_genome --readFilesIn $input_dir/$read1 --readFilesCommand zcat --outFileNamePrefix $output_dir/${prefix}. --runThreadN 10
	echo "Aligning paired-ended reads." 
else 
	echo "Run mode is neither [s]ingle or [p]aired-end." ; exit 1
fi

STAR --genomeDir $star_genome --genomeLoad Remove
echo "FINISH: $(date)" | tee -a $log 
