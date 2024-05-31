#!/usr/bin/env python

import argparse, pandas as pd
from scipy.stats import poisson

# command line parsing info
parser = argparse.ArgumentParser(description='Call PRDM9 alleles')
parser.add_argument('-r', '--reads', type=argparse.FileType('r'), required=True, help='kmers from sequencing reads')
parser.add_argument('-a', '--allele', type=argparse.FileType('r'), required=True, help='kmers from all alleles')
parser.add_argument('-l', '--lam', type=float, default=9, help='lambda for poisson distribution')
parser.add_argument('-e', '--lam_error', type=float, default=2, help='lambda error for sequencing errors')
#parser.add_argument('-o', '--out', type=str, default='out', help='name of output file (defaults to "out")')
args = parser.parse_args()

def count_score(lam=args.lam, lam_error=args.lam_error):
    with args.allele as file:
        allele_kmers = pd.read_csv(args.allele, header=None, names=['Allele', 'kmer', 'AlleleCount'])
    with args.reads as file:
        reads_kmers = pd.read_csv(args.reads, header=None, names=['Allele', 'kmer', 'ReadCount'])

    sim_allele = reads_kmers['Allele'][0]
    reads_kmers = reads_kmers.drop('Allele', axis=1)

    all_alleles = sorted(list(set(allele_kmers['Allele'])))
    all_kmers = sorted(list(set(allele_kmers['kmer'])))
    print("all possible allele kmers:")
    print(all_kmers)
    print("")

    all_kmers_df = pd.DataFrame({'kmer': all_kmers})

    # iterate over alleles
    for a in all_alleles:
        test_kmers = allele_kmers.loc[allele_kmers['Allele']==a]
        print("kmers being tested (allele being tested):")
        print(test_kmers)
        print("")
        combine = pd.concat([all_kmers_df.set_index('kmer'), test_kmers.set_index('kmer')], axis=1)
        combine = combine.drop('Allele', axis=1)

        # add reads kmers
        combine_reads = pd.concat([combine, reads_kmers.set_index('kmer')], axis=1, join_axes=[combine.index]).fillna(0)
        print("dataframe being used:")
        print(combine_reads)
        print("")

        score = 0

        for k in all_kmers:
            n = combine_reads.loc[k, 'ReadCount']
            c = combine_reads.loc[k, 'AlleleCount']
            print("read count for kmer:")
            print(n)
            print("allele count for kmer:")
            print(c)
            print("")
            print("score for kmer:")

            if c==0:
                # could be from a sequencing error
                score += poisson.logpmf(n, lam_error)
                print(poisson.logpmf(n, lam_error))
            else:
                score += poisson.logpmf(n, c*lam)
                print(poisson.logpmf(n, c*lam))
            print("")

        print("final score  for allele:")
        print(sim_allele + ',' + a + ',' + str(score))
        print("")
        print("---")
        print("")

    # outname = args.out + '.csv'

    # out_df = pd.DataFrame(score, index=names).sort_values(0, ascending=False)
    # out_df.to_csv(outname, header=None)

count_score()
