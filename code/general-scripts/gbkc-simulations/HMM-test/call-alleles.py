#!/usr/bin/env python

import sys, argparse, math, numpy as np, pandas as pd
from scipy.special import logsumexp

# command line parsing info
parser = argparse.ArgumentParser(description='Call PRDM9 alleles')
parser.add_argument('-r', '--reads', type=argparse.FileType('r'), required=True, help='sequencing reads to be aligned')
parser.add_argument('-a', '--alleles', type=argparse.FileType('r'), required=True, help='alleles to which reads will be aligned')
parser.add_argument('-g', '--gap', type=float, default=0.000001, help='gap opening probability (defaults to 0.000001)')
parser.add_argument('-e', '--extend', type=float, default=0.000001, help='gap extension probability (defaults to 0.000001)')
parser.add_argument('-t', '--terminate', type=float, default=0.01, help='termination probability (defaults to 0.01)')
parser.add_argument('-m', '--mismatch', type=float, default=0.001, help='mismatch probability (defaults to 0.001)')
parser.add_argument('-o', '--out', type=str, default='out', help='name of output file (defaults to "out")')
args = parser.parse_args()

def reverse_complement(seq):
	seq_upper = seq.upper()
	seq_rc = seq[::-1].replace('A', 't').replace('C', 'g').replace('G', 'c').replace('T', 'a').upper()
	return seq_rc

def call_alleles(reads=args.reads, alleles=args.alleles, g=args.gap, e=args.extend, t=args.terminate, mme=args.mismatch):
	with args.reads as file:
		reads = file.read().splitlines()

	with args.alleles as file:
		alleles_names = pd.read_csv(args.alleles, header=None)

	# prep input files
#	reads = filter(None, reads)
	alleles = list(alleles_names.loc[:][1])
	names = list(alleles_names.loc[:][0])

	# convert transition and emission probabilities to log space
	bm_trans = np.log(1.-g-t)
	xm_trans = np.log(1.-e-t)
	ym_trans = np.log(1.-e)
	mm_trans = np.log(1.-(2.*g)-t)

	g = np.log(g)
	e = np.log(e)
	t = np.log(t)

	match = np.log(1.-mme)
	mismatch = np.log(mme)

	# store full probabilities of alignment
	all_probs = np.zeros((len(reads), len(alleles)))

	# TODO
	# include alignment to reverse complement

	for a in range(len(alleles)):
		for r in range(len(reads)):
	
			# iterate over forward and reverse complement of read
			best_prob = -math.inf
			best_orientation = -1
			for read_orientation, sequence in enumerate([reads[r], reverse_complement(reads[r])]):

				# generate dynamic programming matrices for forward algorithm
				f_B = np.full((len(reads[r])+1, len(alleles[a])+1), -np.inf)
				f_M = f_B.copy()
				f_X = f_B.copy()
				f_Y = f_B.copy()

				# initiate, allowing read to align anywhere along allele
				f_B[0, :] = 0.

				# dynamic programming
				for i in range(1, len(reads[r])+1):
					for j in range(1, len(alleles[a])+1):
						pxy = match if reads[r][i-1] == alleles[a][j-1] else mismatch

						f_M[i, j] = pxy + logsumexp([(bm_trans + f_B[i-1, j-1]), (mm_trans + f_M[i-1, j-1]),  
							(xm_trans + f_X[i-1, j-1]), (ym_trans + f_Y[i-1, j-1])])

						# qxi and qyj (emission probabilities) are both always 1, so don't include to simplify
						f_X[i, j] = logsumexp([(g + f_B[i-1, j]), (g + f_M[i-1, j]), (e + f_X[i-1, j])])

						f_Y[i, j] = logsumexp([(g + f_M[i, j-1]), (e + f_Y[i, j-1])])

				# terminate, taking most probable position in allele of read-end alignment
				f_E = logsumexp([(t + f_B[-1]), (t + f_M[-1]), (t + f_X[-1])], axis=0)

				# check probability
				if max(f_E) > best_prob:
					best_prob = max(f_E)
					best_orientation = read_orientation
				
			all_probs[r, a] = best_prob
			print("Read %d, best orientation %d" % (r, best_orientation), file=sys.stderr)

	all_sums = np.sum(all_probs, 0)

	outname = args.out + '.csv'

	out_df = pd.DataFrame(all_sums, index=names).sort_values(0, ascending=False)
	out_df.to_csv(outname, header=None)

call_alleles()
