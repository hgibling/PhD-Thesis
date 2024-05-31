#!/usr/bin/env python

import argparse, csv, pandas as pd

# command line parsing info
parser = argparse.ArgumentParser(description='Get called alleles')
parser.add_argument('-s', '--scores', type=argparse.FileType('r'), required=True, help='list of scores for tested alleles')
parser.add_argument('-o', '--out', type=str, default='out', help='name of correct called output file (defaults to "out")')
args = parser.parse_args()

def analyze_scores():
	with args.scores as file:
		scores = pd.read_csv(args.scores, header=None, names=['Simulated', 'Tested', 'Score'])

	scores = scores[(scores.Simulated != 'N') & (scores.Tested != 'N')].assign(Correct=scores.Simulated == scores.Tested)
	max_scores = scores[scores.groupby(['Simulated'])['Score'].transform(max) == scores['Score']]
	max_scores = max_scores.assign(Duplicated=max_scores.Simulated.isin(max_scores[max_scores.duplicated('Simulated')].Simulated))

	outname = args.out + '.csv'

	max_scores.to_csv(outname, header=False, index=False)

analyze_scores()
