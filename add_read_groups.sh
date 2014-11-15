#!/bin/bash 

# Usage: 
# arg1	sorted bam folder 


# set global variable: 
picard=/srv/gs1/software/picard-tools/1.111
hg19=/srv/gs1/projects/montgomery/bliu2/shared/ucsc_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa

# read command line arguments: 
input_dir=$1
output_dir=$1 
log=$output_dir/$(basename $0 .sh)

# create array to store all sorted bam file names
cd $input_dir
inputs=(*.sorted.bam)


echo -e $* | tee $log 
echo "START: $(date)" | tee -a $log

for input in ${inputs[*]}; do 

	output=${input/sorted.bam/addrg.bam}
	echo $input
	java -Xmx4g -jar $picard/AddOrReplaceReadGroups.jar INPUT=$input OUTPUT=$output RGID=group1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=sample1 

done 

echo "FINISH: $(date)" | tee -a $log 