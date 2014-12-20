#!/bin/bash

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
bt=/srv/gs1/projects/montgomery/bliu2/bioinformatics_toolbox
brt=/srv/gs1/projects/montgomery/bliu2/brt
vcftools=/srv/gs1/software/vcftools/0.1.12/bin
# load modules:  
module load samtools/1.1
module load tabix/0.2.6
module load java
module load vcftools/0.1.12
# load functions: 
sh $bt/assert.sh 

wd=/srv/gs1/projects/montgomery/bliu2/ancestry/data/Bosh_ASW_mmPCR_2/data/2014-12-04/ASW_genotypes
individual_list=$wd/individual_list.txt # list of individuals. 
position_list=$wd/2popaims_wlrld2M_150.aims.snpinfo.txt # list of aims sites. 
TGP=/srv/gsfs0/projects/montgomery/shared/phase3v5

# run vcftools 
for i in 19; do
	input=$TGP/ALL.chr$i.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz # 1000 genomes phase 3 vcf file
	output_prefix=$wd/ASW.AIMS.chr$i.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes
	vcftools --gzvcf $input --out $output_prefix --recode --keep $individual_list --positions $position_list
done 