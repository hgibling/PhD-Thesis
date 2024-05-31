#!/usr/bin/env python

import sys

def reverse_complement(sequence_file):
	if sys.argv[1] is "-":
		file = sys.stdin
	else:
		file = open(sequence_file, 'r')
	
	for read in file:
		read = read.strip()
		print(read[::-1].replace('A', 't').replace('C', 'g').replace('G', 'c').replace('T', 'a').upper())

reverse_complement(sys.argv[1])