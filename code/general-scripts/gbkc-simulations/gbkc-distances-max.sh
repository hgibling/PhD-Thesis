#!/bin/bash
#$ -cwd

# command use: ./batch-sim-analyze.sh $1 allele list $2 coverage $3 read length $4 error rate $5 number of iterations $6 fragment length $7 kmer size $8 method $9 save directory


# load modules
module load python/2.7.11
source /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/hmm-venv/bin/activate
source ~/scripts/error-handling.sh


# define directories
basedir=/.mounts/labs/awadallalab/private/hgibling
prdm9dir=$basedir/PRDM9-Project/HMM
codedir=$prdm9dir/PRDM9-Allele-Calling/code

mkdir -p $prdm9dir/simulations/$9
simdir=$prdm9dir/simulations/$9
gbkcdir=$basedir/gbkc/src


# make directories if they don't already exist
mkdir -p $simdir/$2/scores
mkdir -p $simdir/$2/pr
mkdir -p $simdir/pr-meta
mkdir -p $simdir/$2/max
mkdir -p $simdir/max-meta


# combine scores and calculate precision and recall overall and for each allele
for i in $(seq 1 $5)
do
	if [ ! -s $simdir/$2/max/all-alleles-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores-max-labeled.csv ]
	then
		$codedir/get-maximum-scores.py -s $simdir/$2/scores/all-alleles-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores.csv -o $simdir/$2/max/all-alleles-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores-max
		awk -v c="$2" -v e="$4" -v k="$7" -v i="$i" -v l="$3" -v f="$6" -v m="$8" -v OFS="," '{print $0, c, e, k, i, l, f, m}' $simdir/$2/max/all-alleles-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores-max.csv > $simdir/$2/max/all-alleles-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores-max-labeled.csv
		if [ ! -s $simdir/$2/max/all-alleles-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores-max-labeled.csv ]; then fail "Max could not be combined for iteration $i"; fi
	fi
	echo "Max $i complete"
done


# combine meta scores and determine how many times alleles called correctly across iterations
if [ ! -s $simdir/max-meta/all-alleles-all-i-flank10k-$3-bp-$6-frag-$2-x-e-$4-F-$7-k-$8-scores-max.csv ]
then
	cat $simdir/$2/max/all-alleles-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-*-F-$7-k-$8-scores-max-labeled.csv > $simdir/max-meta/all-alleles-all-i-flank10k-$3-bp-$6-frag-$2-x-e-$4-F-$7-k-$8-scores-max.csv
	if [ ! -s $simdir/max-meta/all-alleles-all-i-flank10k-$3-bp-$6-frag-$2-x-e-$4-F-$7-k-$8-scores-max.csv ]; then fail "Could not combine meta scores"; fi
fi

echo "Success!"
