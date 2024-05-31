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
graphdir=$prdm9dir/graphs/allele-ref-graphs/updated
alleledir=$prdm9dir/indiv-alleles/flank1k
savedir=$scratchdir/counts-diploid-python-200x/100


lam=$($codedir/calculate-lambda.py 100 51 100 0)
echo "lam is $lam"

g=$1
a1=$(echo $g | cut -d- -f1)
a2=$(echo $g | cut -d- -f2)


for i in $(seq 1 25)
do
    if [ ! -s $savedir/kmers/$1-flank1k-100-bp-250-frag-100-x-e-0-i-$i.kmercounts  ]
    then
        $codedir/get-kmers.py $savedir/reads/$1-flank1k-100-bp-250-frag-100-x-e-0-i-$i.reads 51 F | sort | uniq -c | awk -v g="$1" -v OFS="," '{print g, $2, $1}' > $savedir/kmers/$1-flank1k-100-bp-250-frag-100-x-e-0-i-$i.kmercounts

        if [ ! -s $savedir/kmers/$1-flank1k-100-bp-250-frag-100-x-e-0-i-$i.kmercounts  ]; then fail "kmer counts not generated for iteration $i"; fi
    fi

    if [ ! -s $savedir/scores/$1-flank1k-100-bp-250-frag-100-x-e-0-i-$i.scores.csv ]
    then
        $codedir/count-score-genotypes.py -r $savedir/kmers/$1-flank1k-100-bp-250-frag-100-x-e-0-i-$i.kmercounts -g $scratchdir/genotype-kmer-counts/all-genotypes-51-k.kmercounts -l $lam -e 1 > $savedir/scores/$1-flank1k-100-bp-250-frag-100-x-e-0-i-$i.scores.csv

        if [ ! -s $savedir/scores/$1-flank1k-100-bp-250-frag-100-x-e-0-i-$i.scores.csv ]; then fail "scores not generated for iteration $i"; fi
    fi


    echo "scores calculated for iteration $i"
done

echo "Success!"
