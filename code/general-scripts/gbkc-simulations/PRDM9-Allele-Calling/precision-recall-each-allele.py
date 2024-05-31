#!/usr/bin/env python

# TODO: incorporate into precision-recall.py

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
	scores = scores.assign(Max = scores.Score.isin(scores.groupby(['Simulated'])['Score'].transform(max)), Type=0)
	scores.loc[(scores['Correct']==True) & (scores['Max']==True), 'Type'] = 'TP'
	scores.loc[(scores['Correct']==False) & (scores['Max']==True), 'Type'] = 'FP'
	scores.loc[(scores['Correct']==True) & (scores['Max']==False), 'Type'] = 'FN'
	scores.loc[(scores['Correct']==False) & (scores['Max']==False), 'Type'] = 'TN'
	scores_results = scores[['Simulated', 'Type', 'Score']].groupby(['Simulated', 'Type']).count().unstack().fillna(0)
	scores_results = scores_results.Score

	if 'TP' not in scores_results.columns.values:
		scores_results = scores_results.assign(TP = 0)
	if 'FP' not in scores_results.columns.values:
		scores_results = scores_results.assign(FP = 0)
	if 'FN' not in scores_results.columns.values:
		scores_results = scores_results.assign(FN = 0)
	if 'TN' not in scores_results.columns.values:
		scores_results = scores_results.assign(TN = 0)	

	scores_results = scores_results.assign(Precision = (scores_results.TP / (scores_results.TP + scores_results.FP)),
		Recall = (scores_results.TP / (scores_results.TP + scores_results.FN)))

	outname = args.out + '.csv'
	scores_results[['Precision', 'Recall']].to_csv(outname, header=False)

analyze_scores()
