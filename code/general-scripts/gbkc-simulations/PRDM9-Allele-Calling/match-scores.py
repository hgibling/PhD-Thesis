#!/usr/bin/env python

import argparse, pandas as pd, numpy as np

parser = argparse.ArgumentParser(description='Find kmers from reads with same scores for tested alleles')
parser.add_argument('-a1', '--allele1', type=argparse.FileType('r'), required=True, help='list of kmer scores from allele1')
parser.add_argument('-a2', '--allele2', type=argparse.FileType('r'), required=True, help='list of kmer scores from allele2')
parser.add_argument('-o', '--out', type=str, default='file', help='name for output files')
args = parser.parse_args()


def matching_scores():
	with args.allele1 as file:
		allele1 = pd.read_csv(args.allele1, header=None, names=['Simulated', 'Tested1', 'ID', 'Score1'])
	with args.allele2 as file:
		allele2 = pd.read_csv(args.allele2, header=None, names=['Simulated', 'Tested2', 'ID', 'Score2'])

	allele2 = allele2.drop(['Simulated'], 1)

	both_alleles = pd.merge(allele1, allele2, how='inner', on='ID')
	matches = both_alleles[both_alleles['Score1'] == both_alleles['Score2']]
	mismatches = both_alleles[both_alleles['Score1'] != both_alleles['Score2']]

	match_name = args.out + '-matches.csv'
	mismatch_name = args.out + '-mismatches.csv'

	matches.to_csv(match_name, header=None, index=False)
	mismatches.to_csv(mismatch_name, header=None, index=False)

matching_scores()

