#!/bin/bash

# set global variable: 
gatk=/srv/gs1/software/gatk/gatk-3.3.0/GenomeAnalysisTK.jar
hg19=/srv/gs1/projects/montgomery/bliu2/shared/ucsc_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa

# read command line arguments: 
bam=$1 
echo $bam
vcf=${bam/addrg.bam/vcf}
echo $vcf 

java -Xmx4g -jar $gatk -l INFO -R $hg19 -T UnifiedGenotyper -I $bam -o $vcf --output_mode EMIT_ALL_SITES -L /srv/gs1/projects/montgomery/bliu2/ancestry/data/autosomes/2popaims_wlrld2M_150.aims.snpinfo.bed