#!/usr/bin/env python

import sys, numpy as np

def simulate_reads(file, n, l):
	n = int(n)
	l = int(l)
	f = open(file, 'r')
	for seq in f:
		seq = seq.strip()
		if ',' in seq:
			seq = seq.split(",")[1]
		print(seq)
		seqlength = len(seq)
		nuc = ['A', 'C', 'G', 'T']

		for i in xrange(n):
			pos = np.random.randint(0, (seqlength-l+1))
			print(seq[pos:pos+l])

simulate_reads(sys.argv[1], sys.argv[2], sys.argv[3])
