#!/usr/bin/env python

import argparse, pandas as pd, numpy as np
from scipy.stats import norm
from scipy.stats.mstats import gmean


# command line parsing info
parser = argparse.ArgumentParser(description='Score read kmers based on distances in alleles')
parser.add_argument('-d', '--allele_distances', type=argparse.FileType('r'), required=True, help='list of all allele kmer distances')
parser.add_argument('-k1', '--kmers_pe1', type=argparse.FileType('r'), required=True, help='list of read1 kmer distances')
parser.add_argument('-k2', '--kmers_pe2', type=argparse.FileType('r'), required=True, help='list of read2 kmer distances')
parser.add_argument('-f', '--fragment', type=float, default=250, help='average fragment length of simulated reads')
parser.add_argument('-s', '--standard_dev', type=float, default=50, help='standard deviation of fragments from simulated reads')
parser.add_argument('-p', '--penalty', type=float, default=10000, help='penalty for kmer pairs not observed in allele')
parser.add_argument('-o', '--output', type=str, default='F', help='output score for each read pair for troubleshooting')
# parser.add_argument('-o', '--out', type=str, default='out', help='name of  output file (defaults to "out")')
args = parser.parse_args()


def score_kmer_distances():
	# get input files and prep
	with args.allele_distances as file:
		allele_distances = pd.read_csv(args.allele_distances, header=None, names=['Allele', 'First', 'Second', 'Length'], dtype=str)
	all_alleles = sorted(set(list(allele_distances.Allele)))
	all_kmer_pairs = allele_distances.drop(['Allele', 'Length'], 1).drop_duplicates()

	with args.kmers_pe1 as file:
		kmers_pe1 = pd.read_csv(args.kmers_pe1, header=None, names=['ID1', 'First', 'Position1'], dtype={'ID1': str, 'First': str, 'Position1': int})
	with args.kmers_pe2 as file:
		kmers_pe2 = pd.read_csv(args.kmers_pe2, header=None, names=['ID2', 'Second', 'Position2'], dtype={'ID1': str, 'Second': str, 'Position2': int})

	penalty = args.penalty


	# clean paired end distances and merge together, removing kmers that don't exist in any allele
	kmers_pe1['Allele'], kmers_pe1['ID'] = kmers_pe1['ID1'].str.split(':').str
	kmers_pe1['Allele'], kmers_pe1['Read'] = kmers_pe1['Allele'].str.split('_').str
	kmers_pe1 = kmers_pe1.drop(['ID1', 'Read'], 1)

	kmers_pe2['Allele'], kmers_pe2['ID'] = kmers_pe2['ID2'].str.split(':').str
	kmers_pe2 = kmers_pe2.drop(['ID2', 'Allele'], 1)

	kmers_distances = pd.merge(kmers_pe1, kmers_pe2, how='outer', on='ID')
	sim_allele = kmers_distances['Allele'][0]


	# remove read kmer pairs from flanking sequences by merging with allele kmers
	kmers_distances_filtered = pd.merge(kmers_distances, all_kmer_pairs, how='right', on=['First', 'Second'])


	# for now, filter out all kmers except those at outer ends of reads
	kmers_distances_filtered = kmers_distances_filtered[(kmers_distances_filtered.Position1==0) & (kmers_distances_filtered.Position2==0)]


	# iterate over each allele
	for a in all_alleles:
		test_allele = allele_distances.loc[allele_distances['Allele']==a].drop('Allele', 1)
		test_allele_matrix = test_allele.set_index('Second').pivot(columns='First', values='Length')
		test_allele_dictionary = test_allele_matrix.to_dict('dict')

		score = 0


		# get distance(s) of kmer pair(s) in reads
		for i in xrange(kmers_distances_filtered.shape[0]):
			if kmers_distances_filtered.iloc[i].First in test_allele_dictionary:
				if kmers_distances_filtered.iloc[i].Second in test_allele_dictionary[kmers_distances_filtered.iloc[i].First]:
					all_distances = test_allele_dictionary[kmers_distances_filtered.iloc[i].First][kmers_distances_filtered.iloc[i].Second].split(':')
					distance = map(float, all_distances)
				else:
					distance = penalty
			else:
				distance = penalty

			# for now, take average of scores for all distances from a kmer pair
			pair_score = np.mean(norm.logpdf(distance, loc=args.fragment, scale=args.standard_dev))
			if args.output == 'T':
				print(','.join([sim_allele, a, kmers_distances_filtered.iloc[i].ID, kmers_distances_filtered.iloc[i].First, kmers_distances_filtered.iloc[i].Second, str(pair_score), str(distance)]))
			score += pair_score


		print(','.join([sim_allele, a, str(score)]))


	# outname = args.out + '.csv'
	# kmers_dataframe.to_csv(outname, header=None, index=False)


score_kmer_distances()
