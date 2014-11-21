#!/bin/bash 

# Usage: 
# 
# 
# 
# MODIFY: 
extension=""
output_extension=""

################# SCG settings ################### 
# Job Name 
#$ -N STAR_genome
# 
# Request Large Memory Machine  
# -P large_mem
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
#$ -pe shm 10                                 
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



# DO NOT MODIFY:
software=/srv/gs1/software/
hg19=/srv/gs1/projects/montgomery/bliu2/shared/ucsc_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
gatk=$software/gatk/gatk-3.3.0/GenomeAnalysisTK.jar
picard=$software/picard-tools/1.111
bwa=/srv/gs1/software/bwa/bwa-0.7.7/bin/bwa
star_genome_v14=/srv/gs1/projects/montgomery/tnance/genomes/STAR/hg19_gencode14_overhang99
star_genome_v21=/srv/gs1/projects/montgomery/bliu2/shared/ucsc_hg19/Homo_sapiens/UCSC/hg19/Sequence/STARGenome/hg19_gencode21_overhang75
module load STAR
htseq=/srv/gs1/projects/montgomery/bliu2/tools/HTSeq-0.6.1/scripts
python=/home/bliu2/anaconda/bin/python
gencode14=/srv/gs1/projects/montgomery/shared/annotation/gencode.v14.annotation.gtf
gencode21=/srv/gs1/projects/montgomery/shared/annotation/gencode.v21.annotation.gtf
module load java 

# create array to store all sorted bam file names


# your command here: 
STAR --runMode genomeGenerate --genomeDir $star_genome_v21 --sjdbOverhang 75 --runThreadN 10 --sjdbGTFfile $gencode21 --genomeFastaFiles $hg19

# report done: 
touch $star_genome_v21/hg19_gencode21_overhang75.done  

