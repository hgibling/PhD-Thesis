#!/bin/bash
#$ -cwd

# command use: ./batch-sim-analyze.sh $1 genotype $2 coverage $3 read length $4 error rate $5 number of iterations $6 fragment length $7 upper kmer size $8 method $9 save directory


# load modules
source ~/scripts/error-handling.sh
source /.mounts/labs/simpsonlab/users/hgibling/venvs/gbkc/bin/activate
module load python


# define directories
basedir=/.mounts/labs/awadallalab/private/hgibling
prdm9dir=$basedir/PRDM9-Project/HMM
codedir=$prdm9dir/PRDM9-Allele-Calling/code

#mkdir -p /.mounts/labs/simpsonlab/users/hgibling/simulations/$9
simdir=$prdm9dir/simulations/distances-diploid
gbkcdir=$basedir/gbkc/src


# make directories if they don't already exist
mkdir -p $simdir/$2/scores
mkdir -p $simdir/$2/pr
mkdir -p $simdir/pr-meta


# combine scores and calculate precision and recall overall and for each allele
for i in $(seq 1 $5)
do
    if [ ! -s $simdir/$2/scores/all-genotypes-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores.csv ]
    then
        cat $simdir/$2/scores/*flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores-full.csv > $simdir/$2/scores/all-genotypes-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores.csv
        if [ ! -s $simdir/$2/scores/all-genotypes-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores.csv ]; then fail "Could not combine scores for iteration $i"; fi
    fi
    if [ ! -s $simdir/$2/pr/all-genotypes-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores-pr-labeled.csv ]
    then
        $codedir/precision-recall-k-range.py -s $simdir/$2/scores/all-genotypes-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores.csv -o $simdir/$2/pr/all-genotypes-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores-pr
        awk -v c="$2" -v e="$4" -v i="$i" -v l="$3" -v f="$6" -v m="$8" -v OFS="," '{print $0, c, e, i, l, f, m}' $simdir/$2/pr/all-genotypes-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores-pr.csv > $simdir/$2/pr/all-genotypes-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores-pr-labeled.csv
        if [ ! -s $simdir/$2/pr/all-genotypes-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores-pr-labeled.csv ]; then fail "Scores could not be combined for iteration $i"; fi
    fi
    echo "Precision and recall for iteration $i complete"
done


# combine meta scores and determine how many times alleles called correctly across iterations
if [ ! -s $simdir/pr-meta/all-genotypes-all-i-flank10k-$3-bp-$6-frag-$2-x-e-$4-F-$7-k-$8-scores-pr.csv ]
then
    cat $simdir/$2/pr/all-genotypes-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-*-F-$7-k-$8-scores-pr-labeled.csv > $simdir/pr-meta/all-genotypes-all-i-flank10k-$3-bp-$6-frag-$2-x-e-$4-F-$7-k-$8-scores-pr.csv
    if [ ! -s $simdir/pr-meta/all-genotypes-all-i-flank10k-$3-bp-$6-frag-$2-x-e-$4-F-$7-k-$8-scores-pr.csv ]; then fail "Could not combine meta scores"; fi
fi

echo "Success!"
