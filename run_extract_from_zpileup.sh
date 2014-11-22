#!/bin/bash 

# Usage: 
# 
# SCG SETTINGS: 
# -t 1:<num of jobx>
# 
# CMD ARGS:
# arg1	input dir
# arg2	output dir
# 
# MODIFY: 
extension="zpileup"
output_extension="150aims.zpileup"

################# SCG settings ################### 
# Job Name 
#$ -N extract_from_zpileup
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
# -pe shm 4                                 
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
input_dir=$1
output_dir=$2 
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
brt=/srv/gs1/projects/montgomery/bliu2/brt
bt=/srv/gs1/projects/montgomery/bliu2/bioinformatics_toolbox
aims=/srv/gs1/projects/montgomery/bliu2/ancestry/data/autosomes/2popaims_wlrld2M_150.aims.snpinfo

# create array to store all sorted bam file names
cd $input_dir
inputs=(*.$extension) # put the file extension here. 
# i=$((SGE_TASK_ID-1))
# echo ${inputs[$i]}
# echo $output

for input in ${inputs[*]}; do 
	output=${input/$extension/$output_extension}

	# $python $bt/extract_from_zpileup.py $aims ${inputs[$i]} $output_dir/$output
	$python $bt/extract_from_zpileup.py $aims $input_dir/$input $output_dir/$output
	touch $output_dir/$output.done  
done 
