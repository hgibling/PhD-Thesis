#!/usr/bin/env python

import sys

def calculate_lambda(length, kmer, coverage, error):
	length, kmer, coverage, error = float(length), float(kmer), float(coverage), float(error)
	lam = (length - kmer + 1) * (coverage / length) * ((1 - error) ** kmer)
	print(lam)

calculate_lambda(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
