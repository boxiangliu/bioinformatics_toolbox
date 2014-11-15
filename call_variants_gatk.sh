#!/bin/bash

# set global variable: 
gatk=/srv/gs1/software/gatk/gatk-3.3.0/GenomeAnalysisTK.jar
hg19=/srv/gs1/projects/montgomery/bliu2/shared/ucsc_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa

# read command line arguments: 
bam=$1 
echo $bam
vcf=${bam/sorted.bam/vcf}
echo $vcf 
java -Xmx4g -jar $gatk -T HaplotypeCaller -R $hg19 -I $bam --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 -o $vcf 

