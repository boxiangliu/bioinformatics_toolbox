#!/bin/bash 

# set global variable: 
picard=/srv/gs1/software/picard-tools/1.111
hg19=/srv/gs1/projects/montgomery/bliu2/shared/ucsc_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa

input=$1 
echo $input
output=${input/addrg.bam/addrg.reordered.bam}
echo $output
java -Xmx4g -jar $picard/ReorderSam.jar INPUT=$input OUTPUT=$output REFERENCE=$hg19
