#!/bin/bash
# ARGS:
# arg1	top-level directory where "SampleSheet.csv" is located 
# arg2	directory where fastq files will go. 

#$ -S  /bin/bash
#$ -j n
#$ -V
#$ -R y
#
# set the name of the job
#$ -N bcl2fastq
#
# set the maximum memory usage
#$ -l h_vmem=16G
#
# set the maximum run time
#$ -l h_rt=24:00:00
#
# notify email
#$ -M bliu2@stanford.edu
#
# send mail when job ends or aborts
#$ -m bae
#
# check for errors in the job submission options
#$ -w e
#
#  set the working directory
# -wd /srv/gs1/projects/montgomery/bliu2/sc/src
#
# set the path to the output
# -o /srv/gs1/projects/montgomery/bliu2/sc/log
#
#  set the path for the error file to be saved
# -e /srv/gs1/projects/montgomery/bliu2/sc/log

input_dir=$1
output_dir=$2

/srv/gs1/software/bcl2fastq2/2.14.01/bcl2fastq --runfolder-dir $input_dir --output-dir $output_dir

# e.g. /srv/gs1/software/bcl2fastq2/2.14.01/bcl2fastq --runfolder-dir /srv/gs1/projects/montgomery/bliu2/sc/raw/141002_NS500418_0017_AH02N0AFXX --output-dir /srv/gs1/projects/montgomery/bliu2/sc/data/fastq

