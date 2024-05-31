#!/usr/bin/env python

import argparse, pandas as pd, numpy as np


# command line parsing info
parser = argparse.ArgumentParser(description='Create distance matrix for kmers')
parser.add_argument('-k', '--kmers', type=argparse.FileType('r'), required=True, help='list of kmers with positions')
parser.add_argument('-o', '--out', type=str, default='out', help='name of  output file (defaults to "out")')
parser.add_argument('-p', '--penalty', type=float, default=10000, help='penalty for kmer pairs not observed in allele')
args = parser.parse_args()


def kmer_distances():
	# get input files and prep
	with args.kmers as file:
		kmer_positions = pd.read_csv(args.kmers, header=None, names=['kmer', 'Position'], dtype={'kmer': str, 'Position': int})

	allele_kmers = list(kmer_positions.kmer)
	k = len(allele_kmers[0])

	penalty = args.penalty


	# make symmetrical matrix of kmer distances ("fragment length" of kmer pairs)
	kmers_matrix = pd.DataFrame(np.zeros((len(allele_kmers), len(allele_kmers))), index=allele_kmers, columns=allele_kmers)
	kmers_matrix.iloc[0] = list(kmer_positions.Position)

	for i in xrange(1, kmers_matrix.shape[0]):
		kmers_matrix.iloc[i] = kmers_matrix.iloc[i-1]-1
	kmers_matrix = kmers_matrix+k

	kmers_matrix = pd.DataFrame(np.triu(kmers_matrix), index=allele_kmers, columns=allele_kmers)
	kmers_matrix[kmers_matrix == 0] = penalty


	# get pairs of kmers
	kmers_matrix.index.name = 'First'
	kmers_dataframe = pd.melt(kmers_matrix.reset_index(), id_vars='First', var_name='Second', value_name='Length')

	kmers_dataframe_penalty = kmers_dataframe.loc[kmers_dataframe.Length == penalty].reset_index(drop=True)
	kmers_dataframe = kmers_dataframe.loc[kmers_dataframe.Length < penalty].reset_index(drop=True)


	# merge distances when kmers present more than once
	kmers_dataframe.Length = kmers_dataframe.Length.astype(str)
	kmers_merged = kmers_dataframe.groupby(['First', 'Second']).apply(lambda x: ':'.join(x['Length'])).reset_index()
	kmers_merged.rename(columns={0: 'Length'}, inplace=True)

	# add back penalty distances, but drop if kmer pair already has a distance
	kmers_final = pd.concat([kmers_merged, kmers_dataframe_penalty]).drop_duplicates(['First', 'Second']) 


	kmers_final.to_csv(args.out, header=None, index=False)


kmer_distances()	