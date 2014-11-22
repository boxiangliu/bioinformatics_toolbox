from collections import Counter, defaultdict
import sys
import re
import argparse

'''
FORMAT:
0 - Chromosome
1 - Position
2 - Type (SNV/INDEL/NOVAR)
3 - Depth (total number of reads)
4 - Coverage (total number of reads that pass filters)
5 - Quality filter failures (total number of reads that fail quality filter)
6 - Reference skips (total number of reference skips)
7 - Reference allele
8 - ** Total reference allele count
9 - Reference read positions (- default)
10 - ** Total alternative allele count
11 - [alt allele, 
12 -  alt alle count, 
13 -  positions]...
'''

# Set up the command line argument parser
parser = argparse.ArgumentParser(description='Parse a samtools pileup file.\nOutput: Chromosome, Position, Type, Depth, FilteredCount, SkippedBases, Ref, RefCount, RefPositions, TotalAltCount, Alt, AltCount, AltPositions [Alt, AltCount, AltPositions]+')
parser.add_argument('-q', default=0, type=int, help='Minimum base quality filter.')
args = parser.parse_args()

# Interpret the minimum base quality filter from command line arguments
base_quality_filter = args.q

'''
Extract and reformat indels from the pileup string.
Uses regex to find all indels, excise them delicately, and reformat the pileup string.
'''
def extract_indels(pileup):
	indel_regex = r'[+-][1-9][0-9]*[ATGCNatgcn]+'
	digits = r'[1-9][0-9]*'
	indels = re.findall(indel_regex, pileup)
	for index, indel in enumerate(indels):
		length = get_first_number_in_string(indel)
		start = len(str(length)) + 1
		full_indel = indel[:start+length]
		short_indel = full_indel[0] + full_indel[2:]
		indels[index] = short_indel
		indel_start = pileup.index(full_indel)
		if pileup[indel_start - 3] == '^':
			full_indel = pileup[indel_start - 2: indel_start + len(full_indel)]
		pileup = pileup.replace(full_indel, '', 1)
		indels[index] = re.sub(digits, '', indels[index])

	if re.search(indel_regex, pileup) is not None:
		print "ERROR: " + pileup,'\n',indels
	assert re.search(indel_regex, pileup) is None
	
	return pileup, indels

# Little utility function to excise the first number in a string and return it as an int (e.g. '...+3AGC..-1G.,.' returns 3)
def get_first_number_in_string(string):
	return int(re.findall(r'\d+', string)[0])

'''
Converts a pileup string into a list of individual calls (groups start reads, end reads, individual reads).
Requires that all indels have been previously removed.
Removes ^ to mark starts of reads, leaves the $ for later filtering.
'''
def convert_pileup_to_list(pileup):
	calls = []
	index = 0

	while index < len(pileup):
		base = pileup[index]
		next_base = pileup[index+1] if index < len(pileup) - 1 else ''

		# group together start and ends of reads
		if base == '^':
			if pileup[index+2] == '^':
				print pileup[index:index+6]
			calls.append(pileup[index+2])
			index += 2
		elif next_base == '$':
			calls.append(pileup[index:index+1])
			index += 1
		else:
			calls.append(base) 

		# advance to the next base otherwise
		index += 1

	return calls

# Go through every line of the file
for line_number, line in enumerate(sys.stdin):
	data = line.strip().split('\t')

	# Skip this line if there were no reads in the pileup
	if len(data) == 4:
		continue

	# Auto-detect if the read positions were included in the pileup
	require_read_positions = True if len(data) == 6 else False
	
	# Append faux read index for splat
	if require_read_positions:
		data.append('-')

	# Splat out the data
	chromosome, position, ref, depth, nucleotides, qualities, read_indices = data
	#print data
	
	# Format some of the input
	ref = ref.upper()
	qualities = list(qualities)
	pileup = nucleotides.upper()

	# Format the read positions
	if require_read_positions:
		read_indices = ','.join(list('-' *len(qualities)))
	read_indices = read_indices.split(',')

	# Initialize aggregators
	positions_failed_base_quality_filter = 0
	reference_allele_count = 0
	skipped_bases = 0
	alternate_alleles = Counter()
	reference_positions = []
	allele_postitions = defaultdict(list)

	# Extract and reformat indels
	pileup = re.sub('\^.', '', pileup)
	pileup, indels = extract_indels(pileup)
	indel_alleles = Counter(indels)
	assert sum(indel_alleles.values()) == len(indels)

	# Convert pileup to a list of calls, then zip with qualities and read positions
	pileup = convert_pileup_to_list(pileup)
	calls = zip(pileup, qualities, read_indices)
	calls = [(call, quality, read_index) for call, quality, read_index in calls if '$' not in call] # filter out end reads
	assert any([len(call) == 2 for call, quality, read_index in calls]) == False

	# Iterate over calls, filter out poor quality reads and classify the call
	for call_data in calls:
		call, quality, read_index = call_data
		if not ord(quality) - 33 >= base_quality_filter:
			positions_failed_base_quality_filter += 1
			continue
		else:
			if call in '.,':
				reference_allele_count += 1
				reference_positions.append(read_index)
			elif call in 'ATGCNatgcn':
				alternate_alleles[call] += 1
				allele_postitions[call].append(read_index)
			elif call in '<>*':
				skipped_bases += 1
			else:
				print "-----"
				print line
				print len(pileup)
				raise RuntimeError("We've got a weird nucleotide error on line %i: Base = %s\tQuality = %s" % (line_number, call, quality))

	# Make sure we found the right alternate amount of reference, alternate, failed, and skipped bases
	coverage = reference_allele_count + sum(alternate_alleles.values())
	if int(depth) != reference_allele_count + sum(alternate_alleles.values()) + positions_failed_base_quality_filter + skipped_bases:
		print "-- DEBUGGING --"
		print line
		print "depth: ", int(depth)
		print "sum: ", reference_allele_count + sum(alternate_alleles.values()) + positions_failed_base_quality_filter + skipped_bases
		print "coverage: ", coverage
		print "ref count: ", reference_allele_count
		print "alt count: ", sum(alternate_alleles.values())
		print "bq filter: ", positions_failed_base_quality_filter
		print "skipped  : ", skipped_bases
		print "indels   : ", sum(indel_alleles.values())
		print indels
		print pileup
		print "qualities: ", qualities
	assert (int(depth) == coverage + positions_failed_base_quality_filter + skipped_bases)

	if require_read_positions:
		reference_positions = ['-']
		for allele in alternate_alleles.keys():
			allele_postitions[allele] = ['-']

	if sum(alternate_alleles.values()) == 0 and sum(indel_alleles.values()) == 0:
		output = [chromosome, position, 'NOVAR', depth, coverage, positions_failed_base_quality_filter, skipped_bases, ref, reference_allele_count, ','.join(reference_positions), sum(alternate_alleles.values())]
		output.extend(['N', '0', '-'])
		print '\t'.join([str(s) for s in output])
	
	if sum(alternate_alleles.values()) > 0:
		output = [chromosome, position, 'SNV', depth, coverage, positions_failed_base_quality_filter, skipped_bases, ref, reference_allele_count, ','.join(reference_positions), sum(alternate_alleles.values())]
		for alternate_allele, count in alternate_alleles.most_common():
			output.extend([alternate_allele, count, ','.join(allele_postitions[alternate_allele])])
		print '\t'.join([str(s) for s in output])
	
	if sum(indel_alleles.values()) > 0:
		output = [chromosome, position, 'INDEL', depth, coverage, positions_failed_base_quality_filter, skipped_bases, ref, reference_allele_count, ','.join(reference_positions), sum(indel_alleles.values())]
		for indel, count in indel_alleles.most_common():
			output.extend([indel, count, '-'])
		print '\t'.join([str(s) for s in output])

	# Make sure that stdout flushes
	sys.stdout.flush()
