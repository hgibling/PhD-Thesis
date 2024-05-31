#!/bin/bash
#$ -cwd

# command use: ./gbkc-counts.sh $1 genotype $2 coverage $3 read length $4 error rate $5 number of iterations $6 fragment length $7 kmer size $8 save directory


# load modules
module load python/2.7.11
source /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/hmm-venv/bin/activate
source ~/scripts/error-handling.sh


# define directories
basedir=/.mounts/labs/awadallalab/private/hgibling
prdm9dir=$basedir/PRDM9-Project/HMM
codedir=$prdm9dir/PRDM9-Allele-Calling/code

#mkdir -p $prdm9dir/simulations/$8
simdir=/scratch2/groups/simpsonlab/hgibling/counts-diploid-200x-flank1k
gbkcdir=$basedir/gbkc/src


# make directories if they don't already exist
mkdir -p $simdir/$2/scores
mkdir -p $simdir/$2/pr
mkdir -p $simdir/pr-meta


# combine scores and calculate precision and recall overall and for each allele
for i in $(seq 1 $5)
do
	if [ ! -s $simdir/$2/scores/all-genotypes-flank1k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-scores.csv ]
	then
		cat $simdir/$2/scores/*flank1k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-scores-full.csv > $simdir/$2/scores/all-genotypes-flank1k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-scores.csv
		if [ ! -s $simdir/$2/scores/all-genotypes-flank1k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-scores.csv ]; then fail "Could not combine scores for iteration $i"; fi
	fi
	if [ ! -s $simdir/$2/pr/all-genotypes-flank1k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-scores-pr-labeled.csv ]
	then
		$codedir/precision-recall-k-range.py -s $simdir/$2/scores/all-genotypes-flank1k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-scores.csv -o $simdir/$2/pr/all-genotypes-flank1k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-scores-pr
		awk -v c="$2" -v e="$4" -v i="$i" -v l="$3" -v f="$6" -v OFS="," '{print $0, c, e, i, l, f}' $simdir/$2/pr/all-genotypes-flank1k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-scores-pr.csv > $simdir/$2/pr/all-genotypes-flank1k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-scores-pr-labeled.csv
		if [ ! -s $simdir/$2/pr/all-genotypes-flank1k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-scores-pr-labeled.csv ]; then fail "Scores could not be combined for iteration $i"; fi
	fi
	echo "Precision and recall for iteration $i complete"
done


# combine meta scores
if [ ! -s $simdir/pr-meta/all-genotypes-all-i-flank1k-$3-bp-$6-frag-$2-x-e-$4-F-$7-k-scores-pr.csv ]
then
	cat $simdir/$2/pr/all-genotypes-flank1k-$3-bp-$6-frag-$2-x-e-$4-i-*-F-$7-k-scores-pr-labeled.csv > $simdir/pr-meta/all-genotypes-all-i-flank1k-$3-bp-$6-frag-$2-x-e-$4-F-$7-k-scores-pr.csv
	if [ ! -s $simdir/pr-meta/all-genotypes-all-i-flank1k-$3-bp-$6-frag-$2-x-e-$4-F-$7-k-scores-pr.csv ]; then fail "Could not combine meta scores"; fi
fi

echo "Success!"