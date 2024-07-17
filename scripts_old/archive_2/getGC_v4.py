#!/usr/bin/python3
scriptname = "getGC_v4.py"

# This script calculates the GC content for each sequence in a fasta file and outputs it in a table along
# with the contig it belongs to and the positions along that contig. There are several differences compared
# to getGC_v3.py, which I won't bother to list here. 

import time # required to output the time at which the script was run
from sys import stdout # this import is need to flush python output to the stdout (instead of leaving it
# in the buffer
from sys import argv # this import is needed in order for the script to handle command line arguments
import socket # this import allows us to access the name of the computer that the script is being run on
from Bio import SeqIO
from collections import OrderedDict

def getGC(sequencestring, windowsize = 0, interval = 0):
	if windowsize == 0:
		windowsize = len(sequencestring)
	if interval == 0:
		interval = windowsize
	output = []
	count = 0
	while(1):
		if (count * interval) >= len(sequencestring):
			break
		startingpoint = count * interval
		endingpoint = (count + 1) * interval
		thissection = sequencestring[startingpoint:endingpoint]
		if len(thissection) == 0:
			break
		# We count both GC and AC, because the presence of Ns would bias the results
		GC_numbers = thissection.count('G') + thissection.count('C') + thissection.count('g') + thissection.count('c')
		AT_numbers = thissection.count('A') + thissection.count('T') + thissection.count('a') + thissection.count('t')
		# If there aren't enough non-N values, we report -1
		if (GC_numbers + AT_numbers) < (windowsize / 2):
			GC_content = -1
		else:
			GC_content = GC_numbers / (GC_numbers + AT_numbers)
		output += [[startingpoint, GC_content]]
		count += 1
	return output

if __name__ == '__main__':

	if len(argv) == 2:
		input_filename = argv[1]
		window_size = 0
		interval_size = 0
		output_filename = 'getGC_v4.csv'
	elif len(argv) == 3:
		input_filename = argv[1]
		window_size = argv[2]
		interval_size = 0
		output_filename = 'getGC_v4.csv'
	elif len(argv) == 4:
		input_filename = argv[1]
		window_size = argv[2]
		interval_size = argv[3]
		output_filename = 'getGC_v4.csv'
	elif len(argv) == 5:
		input_filename = argv[1]
		window_size = argv[2]
		interval_size = argv[3]
		output_filename = argv[4]
	else:
		raise Exception("Fail. There should be one to four command line arguments (fasta filename [, window_size [, interval_size [, output_filename]]])")


	print('Running ' + scriptname  + ' at ' + time.strftime("%H:%M") + ' on ' + time.strftime("%d/%m/%Y") + ' using machine ' + socket.gethostname())
	print('\tinput_filename: ' + input_filename)
	print('\twindow_size: ' + window_size)
	print('\tinterval_size: ' + interval_size)
	print('\toutput_filename: ' + output_filename)
	print('\n')
	stdout.flush()

	output_file = open(output_filename, 'w')

	all_sequences = SeqIO.parse(input_filename, 'fasta')
	for seq in all_sequences:
		print('Calculating GC content for contig ' + seq.name)
		this_output = getGC(seq.seq, int(window_size), int(interval_size))
		print('Writing output for ' + seq.name + ' to ' + output_filename)
		for out in this_output:
			output_file.write('\t'.join([seq.name, str(out[0]), str(out[1])]) + '\n')

	output_file.close()

	print('\n\nScript finished running at ' + time.strftime("%H:%M") + ' on ' + time.strftime("%d/%m/%Y") + '\n')
	stdout.flush()

