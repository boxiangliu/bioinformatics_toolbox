# INPUT
# arg1	a directory where sam files are located. It is okay if other files are also in the folder. The script will look for *.sam exntension. 

module load samtools/0.1.19

sam_dir=$1
bam_dir=$1 
echo "INPUT: $sam_dir" > $sam_dir/compress_and_sort_sam.log
echo "START: $(date)"

cd $sam_dir
sam_files=(*.sam)
echo ${#sam_files[*]}
echo ${sam_files[*]}

for sam_file in ${sam_files[*]}; do 
	echo "START: $(date)" >> $sam_dir/compress_and_sort_sam.log
	echo "INPUT: $sam_file" >> $sam_dir/compress_and_sort_sam.log

 	bam_file=${sam_file/.sam/.sorted.bam} # get filename without extension 

 	echo "OUTPUT: $bam_file" >> $sam_dir/compress_and_sort_sam.log

	samtools view -bS $sam_dir/$sam_file | samtools sort -f - $sam_dir/$bam_file
	
	echo -e "FINISH: $(date)\n" >> $sam_dir/compress_and_sort_sam.log

done

echo "FINISH: $(date)"