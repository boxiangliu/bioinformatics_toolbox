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
# -t 	target sites file 
# MODIFY: 
input_extension="bam"
output_extension="vcf"

################# SCG settings ################### 
# Job Name 
#$ -N unifiedGenotyper
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
while getopts ":i:o:t:" opt; do 
	case $opt in 
		i) input_dir=$OPTARG
		[ $parrot ] & echo "INPUT: $input_dir"
		;;
		o) output_dir=$OPTARG
		[ $parrot ] & echo "OUTPUT: $output_dir"
		;;
		t) target_sites_file=$OPTARG
		[ $parrot ] & echo "TARGET SITES: $target_sites_file"
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
inputs=(*$input_extension) # put the file extension here. 
i=$((SGE_TASK_ID-1))
output=${inputs[$i]/$input_extension/$output_extension}

start=$(date)

java -Xmx4g -jar $gatk -R $hg19 -T UnifiedGenotyper -I $input_dir/${inputs[$i]} -o $output_dir/$output --output_mode EMIT_ALL_SITES -L $target_sites_file

finish=$(date)
touch $output_dir/$output.done; echo -e "START: $start\nFINISH: $finish" > $output_dir/$output.done

################ qlogin script, for debugging ################

#!/bin/bash

# set global variable: 
# gatk=/srv/gs1/software/gatk/gatk-3.3.0/GenomeAnalysisTK.jar
# hg19=/srv/gs1/projects/montgomery/bliu2/shared/ucsc_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa

# # read command line arguments: 
# bam=$1 
# echo $bam
# vcf=${bam/sorted.bam/vcf}
# echo $vcf 

# java -Xmx4g -jar $gatk -l INFO -R $hg19 -T UnifiedGenotyper -I $bam -o $vcf --output_mode EMIT_ALL_SITES -L /srv/gs1/projects/montgomery/bliu2/ancestry/data/autosomes/2popaims_wlrld2M_150.aims.snpinfo.bed