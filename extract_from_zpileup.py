# -*- coding: utf-8 -*-
"""
Created on Thu Nov 20 21:40:53 2014

@author: boshliu

Usage:
arg1 file1 with sites of interest
arg2 zpileup file 
arg3 zpileup file with only sites of interest

Modify the format that describes the position of chromosome and position in the target_sites file. 
"""
target_sites_format = {'chr_field': 0, 'pos_field': 1}

import sys 
target_sites_filename = sys.argv[1]
zpileup_filename = sys.argv[2]
output_filename = sys.argv[3]
target_sites = open(target_sites_filename, 'r')
zpileup = open(zpileup_filename, 'r')
output = open(output_filename, 'w')
lines_to_find = []

zpileup_locations = {}
for line in zpileup:
	split_line = line.strip().split()
	chr = split_line[0][3:]
	pos = split_line[1]
	zpileup_locations[(chr,pos)] = line 


for line in target_sites:
	if 'position' in line.strip(): continue 
	split_line = line.strip().split()
	chr = split_line[target_sites_format['chr_field']]
	if 'chr' in chr: chr = chr[3:] 
	pos = split_line[target_sites_format['pos_field']]
	if (chr, pos) in zpileup_locations:
		output.write(zpileup_locations[(chr,pos)])
	else: 
		output.write("chr%s\t%s\n"%(chr,pos))

output.close()

# for line in target_sites:
#     split_line = line.strip().split()
#     print(split_line)
#     lines_to_find.append((split_line[1], split_line[2]))


# zpileup = zpileup.readlines()
# for line in zpileup:
#     split_line = line.strip().split()
#     chr = split_line[0][3:]
#     pos = split_line[1]
#     print((chr, pos))
#     if (chr, pos) in lines_to_find:

#         output.write(line)

# output.close() 



