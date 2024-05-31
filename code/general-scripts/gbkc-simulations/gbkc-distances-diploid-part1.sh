#!/bin/bash
#$ -cwd

# command use: ./gbkc-counts.sh
# $1 allele
# $2 coverage
# $3 read length
# $4 error rate
# $5 number of iterations
# $6 fragment length
# $7 upper kmer size
# $8  method
# $9 save directory


# load modules
source ~/scripts/error-handling.sh


# define directories
basedir=/.mounts/labs/awadallalab/private/hgibling
prdm9dir=$basedir/PRDM9-Project/HMM
codedir=$prdm9dir/PRDM9-Allele-Calling/code
gbkcdir=$basedir/gbkc/src
genotypedir=$prdm9dir/individual-genotypes/flank10k
alleledir=$prdm9dir/individual-haplotypes/flank10k
simdir=$prdm9dir/simulations/distances-diploid


# make directories if they don't already exist
mkdir -p $simdir/$2/scores
mkdir -p $simdir/$2/pr
mkdir -p $simdir/pr-meta


# iterate over number of simulations
g=$1
geno=$(echo $g | sed 's/-/\//')

for i in $(seq 1 $5)
do
    # calculate scores
    if [ ! -s $simdir/$2/scores/$g-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores-full.csv ]
    then
        $gbkcdir/gbkc distance -a $prdm9dir/all-alleles-flank10k.fa -1 $simdir/$2/reads/$g-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i.r1.fq -2 $simdir/$2/reads/$g-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i.r2.fq -K $7 -l $3 -e $4 -c $2 -f $6 -s 50 -m $8 -d -o $simdir/$2/scores/$g-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores.csv
	if [ $? -ne 0 ]; then fail "gbkc aborted during iteration $i"; fi
        awk -v g="$geno" -v FS="," -v OFS="," '{print $1, g, $2, $3}' $simdir/$2/scores/$g-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores.csv > $simdir/$2/scores/$g-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores-full.csv
        if [ ! -s $simdir/$2/scores/$g-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores-full.csv ]; then fail "No scores for genotype $g iteration $i"; fi
	rm $simdir/$2/scores/$g-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores.csv
    fi
    echo "Scoring for iteration $i complete"
done

echo "Success!"
