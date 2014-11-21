# -*- coding: utf-8 -*-
"""
Created on Thu Nov 20 21:40:53 2014

@author: boshliu
"""

import sys 
source_filename = sys.argv[1]
target_filename = sys.argv[2]

source = open(source_filename, 'r')
target = open(target_filename, 'r')
lines_to_find = []
for line in source:
    split_line = line.strip().split()
    print(split_line)
    lines_to_find.append((split_line[1], split_line[2]))
output = []
target = target.readlines()
for line in target:
    split_line = line.strip().split()
    chr = split_line[0][3:]
    pos = split_line[1]
    print((chr, pos))
    if (chr, pos) in lines_to_find:
        output.append(line)

sys.stdout.write(output)


