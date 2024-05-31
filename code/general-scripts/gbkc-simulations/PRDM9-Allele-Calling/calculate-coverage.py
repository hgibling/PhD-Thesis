#!/usr/bin/env python

import sys, csv

def calculate_coverage(file, coverage, read_length, paired_end):

	with open(file, 'rb') as csvfile:
		reader = csv.reader(csvfile)
		for row in reader:
			sequence = row[1]
			length = float(len(sequence))
			coverage = float(coverage)
			read_length = float(read_length)
			count = (coverage * length) / read_length
			
			if paired_end == 'T':
				print(int(round(count / 2)))
			else:
				print(int(round(count)))


calculate_coverage(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
