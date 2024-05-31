#!/bin/bash
#$ -cwd

# command use: ./batch-sim-analyze.sh 
# $1 NA 
# $2 coverage 
# $3 read length 
# $4 error rate 
# $5 number of iterations 
# $6 fragment length 
# $7 NA 
# $8 method 
# $9 NA


# load modules
module load python/2.7.11
source /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/hmm-venv/bin/activate
source ~/scripts/error-handling.sh


# define directories
basedir=/.mounts/labs/awadallalab/private/hgibling
prdm9dir=$basedir/PRDM9-Project/HMM
codedir=$prdm9dir/PRDM9-Allele-Calling/code

countdir=$prdm9dir/simulations/counts-haploid
distdir=$prdm9dir/simulations/distances-haploid
simdir=$prdm9dir/simulations/combine-counts-distances-haploid

# make directories if they don't already exist
mkdir -p $simdir/$2/scores
mkdir -p $simdir/$2/pr
mkdir -p $simdir/pr-meta

# combine scores and calculate precision and recall overall and for each allele
for i in $(seq 1 $5)
do
    if [ ! -s $simdir/$2/scores/all-alleles-all-k-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores-combined.csv ]
    then
        $codedir/combine-haploid-counts-distances-scores.py -d $distdir/$2/scores-all-k/all-alleles-all-k-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores.csv -c $countdir/$2/scores-all-k/all-alleles-all-k-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-scores.csv -o $simdir/$2/scores/all-alleles-all-k-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores-combined
        if [ ! -s $simdir/$2/scores/all-alleles-all-k-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores-combined.csv ]; then fail "Counts and distances scores could not be combined for iteration $i"; fi
    fi

    if [ ! -s $simdir/$2/pr/all-alleles-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores-pr-labeled.csv ]
    then
        $codedir/precision-recall-k-range.py -s $simdir/$2/scores/all-alleles-all-k-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores-combined.csv -o $simdir/$2/pr/all-alleles-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores-pr
        awk -v c="$2" -v e="$4" -v i="$i" -v l="$3" -v f="$6" -v m="$8" -v OFS="," '{print $0, c, e, i, l, f, m}' $simdir/$2/pr/all-alleles-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores-pr.csv > $simdir/$2/pr/all-alleles-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores-pr-labeled.csv
        if [ ! -s $simdir/$2/pr/all-alleles-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores-pr-labeled.csv ]; then fail "Precision and recall could not be completed for iteration $i"; fi
    fi
    echo "Precision and recall for iteration $i complete"
done


# combine meta scores and determine how many times alleles called correctly across iterations
if [ ! -s $simdir/pr-meta/all-alleles-all-i-flank10k-$3-bp-$6-frag-$2-x-e-$4-F-$7-k-$8-scores-pr.csv ]
then
    cat $simdir/$2/pr/all-alleles-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-*-F-$7-k-$8-scores-pr-labeled.csv > $simdir/pr-meta/all-alleles-all-i-flank10k-$3-bp-$6-frag-$2-x-e-$4-F-$7-k-$8-scores-pr.csv
    if [ ! -s $simdir/pr-meta/all-alleles-all-i-flank10k-$3-bp-$6-frag-$2-x-e-$4-F-$7-k-$8-scores-pr.csv ]; then fail "Could not combine meta scores"; fi
fi

echo "Success!"