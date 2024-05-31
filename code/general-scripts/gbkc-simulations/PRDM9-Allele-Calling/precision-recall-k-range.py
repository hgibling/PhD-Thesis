#!/usr/bin/env python

import argparse, csv, pandas as pd

# command line parsing info
parser = argparse.ArgumentParser(description='Get called alleles')
parser.add_argument('-s', '--scores', type=argparse.FileType('r'), required=True, help='list of scores for tested alleles')
parser.add_argument('-o', '--out', type=str, default='out', help='name of correct called output file (defaults to "out")')
args = parser.parse_args()

def analyze_scores():
	with args.scores as file:
		scores = pd.read_csv(args.scores, header=None, names=['k', 'Simulated', 'Tested', 'Score'])

	scores = scores[(scores.Simulated != 'N') & (scores.Tested != 'N')].assign(Correct=scores.Simulated == scores.Tested)
	scores = scores.assign(Max = scores.Score.isin(scores.groupby(['k', 'Simulated'])['Score'].transform(max)), Type=0)
	scores.loc[(scores['Correct']==True) & (scores['Max']==True), 'Type'] = 'TP'
	scores.loc[(scores['Correct']==False) & (scores['Max']==True), 'Type'] = 'FP'
	scores.loc[(scores['Correct']==True) & (scores['Max']==False), 'Type'] = 'FN'
	scores.loc[(scores['Correct']==False) & (scores['Max']==False), 'Type'] = 'TN'
	scores_results = scores[['k', 'Simulated', 'Type', 'Score']].groupby(['k', 'Simulated', 'Type']).count().unstack().fillna(0)
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
		Recall = (scores_results.TP / (scores_results.TP + scores_results.FN))).fillna(0)

	scores_results = scores_results.assign(F1Score = 2 * ((scores_results.Precision * scores_results.Recall) / (scores_results.Precision + scores_results.Recall))).fillna(0)

	scores_results.reset_index(level=['k', 'Simulated'], inplace=True)

	scores_averages = scores_results.groupby(['k']).mean()
	scores_averages['Simulated'] = "Average"
	scores_averages.reset_index(level=['k'], inplace=True)

	scores_summary = scores_results.append(scores_averages)
	scores_summary = scores_summary[['k', 'Simulated', 'Precision', 'Recall', 'F1Score']].fillna(0)

	outname = args.out + '.csv'

	scores_summary.to_csv(outname, header=False, index=False)

analyze_scores()
