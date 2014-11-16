#!/bin/sh 
#
# set the name of the job 
#$ -N STAR 2-pass
# 
# tell the scheduler this takes a lot of memory 
#$ -P large_mem
# 
# set the maximum memory usage 
#$ -l h_vmem=40G 
#
# set the maximum run time 
#$ -l h_rt=96:00:00  
# 
# check for errors in the job submission options
#$ -w e
#
# email when job ends or aborts
#$ -m ea
#
# specify email address
#$ -M jollier.liu@gmail.com
# 
# run on multiple threads                     
#$ -pe shm 4                                 
#
# run job in current working directory      
#$ -cwd                                    
#
# run job in specified working directory
# -wd /srv/gs1/projects/montgomery/haglund/genome/

# set the output file
# -o  whatever.log
#
# merge standard error stream into the standard output stream                                            
# -j y

# Usage: 
# arg1	fastq directory
# arg2	alignment directory 

set -u 
set -e 

# read command line arguments: 
input_dir=$1
cd $input_dir
output_dir=$2
log=$output_dir/$(basename $0 .sh).log

# modify the suffix for R1 and R2:
fastq_R1=(*_R1_001.fastq.gz)
fastq_R2=(*_R2_001.fastq.gz)

# DO NOT MODIFY:
software=/srv/gs1/software/
hg19=/srv/gs1/projects/montgomery/bliu2/shared/ucsc_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
gatk=$software/gatk/gatk-3.3.0/GenomeAnalysisTK.jar
picard=$software/picard-tools/1.111
bwa=/srv/gs1/software/bwa/bwa-0.7.7/bin/bwa
star_genome=/srv/gs1/projects/montgomery/tnance/genomes/STAR/hg19_gencode14_overhang99
module load STAR/2.3.0e

#### main: 

# check the number of forward and reverse files match:
num_R1=${#fastq_R1[*]}
num_R2=${#fastq_R2[*]}
[ $num_R1 -eq $num_R2 ] || (echo "Number of forward read files and reverse read files not equal"; exit 1)

echo -e "COMMAND: $0 $*" | tee $log
echo "START: $(date)" | tee -a $log

# load genome into shared memory:
STAR --genomeDir $star_genome --genomeLoad LoadAndExit; echo "FINISH LOADING GENOME $(date)" | tee -a $log 

cd $input_dir
for read1 in ${fastq_R1[*]}; do 
	read2=${read1/R1/R2}
	prefix=$(basename $read1 _R1_001.fastq.gz)
	echo -e "SAMPLES: $read1\t$read2" 
	echo "OUTPUT: $output_dir/${prefix}"
	STAR --genomeDir $star_genome --readFilesIn $read1 $read2 --readFilesCommand zcat --outFileNamePrefix $output_dir/${prefix}. --runThreadN 10
done

STAR --genomeDir $star_genome --genomeLoad Remove
echo "FINISH: $(date)" | tee -a $log 
