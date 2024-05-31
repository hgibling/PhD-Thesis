#!/usr/bin/env python

import sys

def get_kmers(sequence_file, k, pos):
	file = open(sequence_file, 'r')		
	k = int(k)
	for read in file:
		read = read.strip()
		if ',' in read:
			read = read.split(",")[1]
		for i in range(len(read)-k+1):
			kmer = read[i:i+k]
			if pos=='T':
				print(",".join([kmer, str(i)]))
			else:
				print(kmer)

get_kmers(sys.argv[1], sys.argv[2], sys.argv[3])
