#!/bin/bash
#$ -cwd

# command use: ./gbkc-counts.sh $1 genotype $2 coverage $3 read length $4 error rate $5 number of iterations $6 fragment length $7 upper kmer size $8 save directory


# load modules
module load python/2.7.11
source /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/hmm-venv/bin/activate
source ~/scripts/error-handling.sh


# define directories
basedir=/.mounts/labs/awadallalab/private/hgibling
prdm9dir=$basedir/PRDM9-Project/HMM
codedir=$prdm9dir/PRDM9-Allele-Calling/code

#mkdir -p $prdm9dir/simulations/$8
simdir=/scratch2/groups/simpsonlab/hgibling/counts-diploid-200x
gbkcdir=$basedir/gbkc/src


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
	if [ ! -s $simdir/$2/scores/$g-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-scores-full.csv ]
	then
		$gbkcdir/gbkc count -a $prdm9dir/all-alleles-flank10k.fa -f $prdm9dir/flank10k.fa -1 $simdir/$2/reads/$g-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i.r1.fq -2 $simdir/$2/reads/$g-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i.r2.fq -K $7 -l $3 -e $4 -c $2 -d -t 12 -o $simdir/$2/scores/$g-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-scores.csv
		awk -v g="$geno" -v FS="," -v OFS="," '{print $1, g, $2, $3}' $simdir/$2/scores/$g-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-scores.csv > $simdir/$2/scores/$g-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-scores-full.csv
		if [ ! -s $simdir/$2/scores/$g-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-scores-full.csv ]; then fail "No scores for allele $g iteration $i"; fi
	fi
	echo "Scoring for iteration $i complete"
done

echo "Success!"