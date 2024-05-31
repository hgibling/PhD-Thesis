#!/usr/bin/env python

import argparse, csv, pandas as pd

# command line parsing info
parser = argparse.ArgumentParser(description='Get called alleles')
parser.add_argument('-s', '--scores', type=argparse.FileType('r'), required=True, help='list of scores for tested alleles')
parser.add_argument('-o', '--out', type=str, default='out', help='name of correct called output file (defaults to "out")')
#parser.add_argument('-l', '--out', type=str, default='list', help='name of called list output file (defaults to "list")')
args = parser.parse_args()

def analyze_scores():
	with args.scores as file:
		scores = pd.read_csv(args.scores, header=None, names=['Simulated', 'Tested', 'Score'])

	scores = scores[(scores.Simulated != 'N') & (scores.Tested != 'N')].assign(Correct=scores.Simulated == scores.Tested)
	max_scores = scores[scores.groupby(['Simulated'])['Score'].transform(max) == scores['Score']]
	max_scores = max_scores.assign(Duplicated=max_scores.Simulated.isin(max_scores[max_scores.duplicated('Simulated')].Simulated))
	max_duplicated = max_scores[max_scores.Duplicated == True].groupby(['Simulated', 'Correct']).count().reset_index().groupby('Simulated').count().reset_index()

	correct = len(max_scores[(max_scores.Correct == True) & (max_scores.Duplicated == False)])
	incorrect = len(max_scores[(max_scores.Correct == False) & (max_scores.Duplicated == False)]) + len(max_duplicated[max_duplicated.Tested == 1])
	tied = len(max_duplicated[max_duplicated.Tested == 2])

	results = [correct, incorrect, tied]

	outname = args.out + '.csv'

	with open(outname, 'wb') as myfile:
		wr = csv.writer(myfile)
		wr.writerow(results)

analyze_scores()
