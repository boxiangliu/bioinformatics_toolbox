#!/bin/bash
# 
# Usage: 
# This script takes genome coordinates in primer3 output and generate a fastq file containing the sequences specified by the coordinates. 
# 
# CMD ARGS:
# -i 	primer3 output file as [i]nput
# -o 	output file 
# 
# 
# Purpose: 
# To generate a reference genome for only target sequences in an mmPCR library. 
# 
# Overview:
# create an empty output file 
# 
# while coordinate in coordinates:
#	read in coordinates of targeted regions  
# 	use samtools faidx to retrieve the targeted sequence 
#	append the coordinates + \n (format: chr#:start-end)
#	append the targeted sequence + \n
# end while 

# Paths: 
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
module load samtools/1.1
sh $bt/assert.sh 

#
# read command line args:
# 
parrot=true 
while getopts ":i:o:" opt; do 
	case $opt in 
		i) input=$OPTARG
		[ $parrot ] & echo "INPUT: $input"
		;;
		o) output=$OPTARG
		[ $parrot ] & echo "OUTPUT: $output"
		;;
		\?) echo "Invalid option -$OPTARG"; exit 1
		;;
	esac
done 
#
# create an empty output file: 
#
if [ ! -f $output ]; then
	echo "Creating $output"
	touch $output
fi 
#
# empty the content of output file. 
# 
if [ -s $output ]; then
	echo -n "" > $output 
fi 
#
# read in coordinates of target regions: 
# 
while read line; do
	[ $parrot ] & echo "Current Line: $line"
	line=($line) # convert string into array
	rs_id=${line[0]}
	coordinate=${line[5]} # format: chr#:start-end
	length=${line[6]}
	#
	# append sequence tag to output file:
	# 
	echo ">$coordinate:$rs_id:$length" >> $output
	#
	# append sequnce to output file:  
	# 
	samtools faidx $hg19 $coordinate | grep -v ">" >> $output
done < $input



