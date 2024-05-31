#!/bin/bash
#$ -cwd

# load modules
module load python/2.7.11
source /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/hmm-venv/bin/activate
source ~/scripts/error-handling.sh

# define directories
basedir=/.mounts/labs/awadallalab/private/hgibling
prdm9dir=$basedir/PRDM9-Project/HMM
codedir=$prdm9dir/PRDM9-Allele-Calling/code
scratchdir=/scratch2/groups/simpsonlab/hgibling
savedir=$scratchdir/counts-diploid-python-200x/100

# combine scores and calculate precision and recall overall and for each allele
for i in $(seq 1 25)
do
    if [ ! -s $savedir/scores/all-genotypes-flank1k-100-bp-250-frag-100-x-e-0-i-$i.scores.csv ]
    then
        cat $savedir/scores/*flank1k-100-bp-250-frag-100-x-e-0-i-$i.scores.csv > $savedir/scores/all-genotypes-flank1k-100-bp-250-frag-100-x-e-0-i-$i.scores.temp.csv
        awk '{print "51," $0}' $savedir/scores/all-genotypes-flank1k-100-bp-250-frag-100-x-e-0-i-$i.scores.temp.csv > $savedir/scores/all-genotypes-flank1k-100-bp-250-frag-100-x-e-0-i-$i.scores.csv
        rm $savedir/scores/all-genotypes-flank1k-100-bp-250-frag-100-x-e-0-i-$i.scores.temp.csv
        if [ ! -s $savedir/scores/all-genotypes-flank1k-100-bp-250-frag-100-x-e-0-i-$i.scores.csv ]; then fail "Could not combine scores for iteration $i"; fi
    fi
    echo "Scores combined for iteration $i"
    if [ ! -s $savedir/pr/all-genotypes-flank1k-100-bp-250-frag-100-x-e-0-i-$i-scores-pr-labeled.csv ]
    then
        $codedir/precision-recall-k-range.py -s $savedir/scores/all-genotypes-flank1k-100-bp-250-frag-100-x-e-0-i-$i.scores.csv -o $savedir/pr/all-genotypes-flank1k-100-bp-250-frag-100-x-e-0-i-$i-scores-pr
        awk -v c="100" -v e="0" -v i="$i" -v l="100" -v f="250" -v OFS="," '{print $0, c, e, i, l, f}' $savedir/pr/all-genotypes-flank1k-100-bp-250-frag-100-x-e-0-i-$i-scores-pr.csv > $savedir/pr/all-genotypes-flank1k-100-bp-250-frag-100-x-e-0-i-$i-scores-pr-labeled.csv
        if [ ! -s $savedir/pr/all-genotypes-flank1k-100-bp-250-frag-100-x-e-0-i-$i-scores-pr-labeled.csv ]; then fail "Scores could not be combined for iteration $i"; fi
    fi
    echo "Precision and recall for iteration $i complete"
done


# combine meta scores
if [ ! -s $savedir/pr-meta/all-genotypes-all-i-flank1k-100-bp-250-frag-100-x-e-0-scores-pr.csv ]
then
    cat $savedir/pr/all-genotypes-flank1k-100-bp-250-frag-100-x-e-0-i-*-scores-pr-labeled.csv > $savedir/pr-meta/all-genotypes-all-i-flank1k-100-bp-250-frag-100-x-e-0-scores-pr.csv
    if [ ! -s $savedir/pr-meta/all-genotypes-all-i-flank1k-100-bp-250-frag-100-x-e-0-scores-pr.csv ]; then fail "Could not combine meta scores"; fi
fi

echo "Success!"
