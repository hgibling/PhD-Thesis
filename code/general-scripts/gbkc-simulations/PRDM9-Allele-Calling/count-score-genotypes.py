#!/usr/bin/env python

import argparse, pandas as pd
from scipy.stats import poisson

# command line parsing info
parser = argparse.ArgumentParser(description='Call PRDM9 genotypes')
parser.add_argument('-r', '--reads', type=argparse.FileType('r'), required=True, help='kmers from sequencing reads')
parser.add_argument('-g', '--genotype', type=argparse.FileType('r'), required=True, help='kmers from all genotypes')
parser.add_argument('-l', '--lam', type=float, default=9, help='lambda for poisson distribution')
parser.add_argument('-e', '--lam_error', type=float, default=2, help='lambda error for sequencing errors')
#parser.add_argument('-o', '--out', type=str, default='out', help='name of output file (defaults to "out")')
args = parser.parse_args()

def count_score(lam=args.lam, lam_error=args.lam_error):
	with args.genotype as file:
		genotype_kmers = pd.read_csv(args.genotype, header=None, names=['Genotype','kmer', 'GenotypeCount'])
	with args.reads as file:
		reads_kmers = pd.read_csv(args.reads, header=None, names=['Genotype', 'kmer', 'ReadCount'])

	sim_genotype = reads_kmers['Genotype'][0]
	reads_kmers = reads_kmers.drop('Genotype', axis=1)

	all_genotypes = sorted(list(set(genotype_kmers['Genotype'])))
	all_kmers = sorted(list(set(genotype_kmers['kmer'])))

	all_kmers_df = pd.DataFrame({'kmer': all_kmers})

	# iterate over genotypes
	for g in all_genotypes:
		test_kmers = genotype_kmers.loc[genotype_kmers['Genotype']==g]
		combine = pd.concat([all_kmers_df.set_index('kmer'), test_kmers.set_index('kmer')], axis=1)
		combine = combine.drop('Genotype', axis=1)

		# add reads kmers
		combine_reads = pd.concat([combine, reads_kmers.set_index('kmer')], axis=1, join_axes=[combine.index]).fillna(0)

		score = 0

		for k in all_kmers:
			n = combine_reads.loc[k, 'ReadCount']
			c = combine_reads.loc[k, 'GenotypeCount']

			if c==0:
				# could be from a sequencing error
				score += poisson.logpmf(n, lam_error)
			else:
				score += poisson.logpmf(n, c*lam)

		print(sim_genotype + ',' + g + ',' + str(score))

	# outname = args.out + '.csv'

	# out_df = pd.DataFrame(score, index=names).sort_values(0, ascending=False)
	# out_df.to_csv(outname, header=None)

count_score()