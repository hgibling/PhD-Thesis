#!/usr/bin/env python

import argparse, pandas as pd, numpy as np
from scipy.stats import norm
from scipy.stats.mstats import gmean


# command line parsing info
parser = argparse.ArgumentParser(description='Score read kmers based on distances in alleles')
parser.add_argument('-d', '--distance_scores', type=argparse.FileType('r'), required=True, help='list of distance scores')
parser.add_argument('-c', '--count_scores', type=argparse.FileType('r'), required=True, help='list of kmer count scores')
parser.add_argument('-o', '--out', type=str, default='out', help='name of  output file (defaults to "out")')
args = parser.parse_args()


def combine_scores():
	# get input files
	with args.distance_scores as file:
		distance_scores = pd.read_csv(args.distance_scores, header=None, names=['Simulated', 'Tested', 'DistanceScore'])
	with args.count_scores as file:
		count_scores = pd.read_csv(args.count_scores, header=None, names=['Simulated', 'Tested', 'CountScore'])

	# combine scores
	combine = pd.merge(distance_scores, count_scores, how='outer', on=['Simulated', 'Tested'])
	combine = combine[combine.Simulated != "N"][combine.Tested != "N"]
	combine['CombinedScore'] = combine.DistanceScore + combine.CountScore


	outname = args.out + '.csv'
	combine.to_csv(outname, header=None, index=False)


combine_scores()
