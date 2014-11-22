# -*- coding: utf-8 -*-
"""
Created on Thu Nov 20 21:40:53 2014

@author: boshliu

Usage:
arg1 file1 with sites of interest
arg2 file2  with all sites
arg3 output with only sites of interest from file2
"""

import sys 
source_filename = sys.argv[1]
target_filename = sys.argv[2]
output_filename = sys.argv[3]
source = open(source_filename, 'r')
target = open(target_filename, 'r')
output = open(output_filename, 'w')
lines_to_find = []

target_locations = {}
for line in target:
	split_line = line.strip().split()
	chr = split_line[0][3:]
	pos = split_line[1]
	target_locations[(chr,pos)] = line 


for line in source:
	if 'position' in line.strip(): continue 
	split_line = line.strip().split()
	chr = split_line[1]
	pos = split_line[2]
	if (chr, pos) in target_locations:
		output.write(target_locations[(chr,pos)])
	else: 
		output.write("chr%s\t%s\n"%(chr,pos))

output.close()

# for line in source:
#     split_line = line.strip().split()
#     print(split_line)
#     lines_to_find.append((split_line[1], split_line[2]))


# target = target.readlines()
# for line in target:
#     split_line = line.strip().split()
#     chr = split_line[0][3:]
#     pos = split_line[1]
#     print((chr, pos))
#     if (chr, pos) in lines_to_find:

#         output.write(line)

# output.close() 



