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


# iterate over list of alleles
while read a
do
	for i in $(seq 1 $5)
	do
		# calculate scores
		if [ ! -s $simdir/$2/scores/$a-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores-full.csv ]
		then
			$gbkcdir/gbkc distance -a $prdm9dir/all-alleles-flank10k.fa -1 $simdir/$2/reads/$a-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i.r1.fq -2 $simdir/$2/reads/$a-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i.r2.fq -k $7 -l $3 -e $4 -c $2 -f $6 -s 50 -m $8 -d -o $simdir/$2/scores/$a-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores.csv
			awk -v a="$a" '{print a "," $0}' $simdir/$2/scores/$a-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores.csv > $simdir/$2/scores/$a-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores-full.csv
			if [ ! -s $simdir/$2/scores/$a-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores-full.csv ]; then fail "No scores for allele $a iteration $i"; fi
		fi
	done
	echo "Scores calculated for allele $a"
done < $1


# combine scores and calculate precision and recall overall and for each allele
for i in $(seq 1 $5)
do
	if [ ! -s $simdir/$2/scores/all-alleles-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores.csv ]
	then
		cat $simdir/$2/scores/*flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores-full.csv > $simdir/$2/scores/all-alleles-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores.csv
		if [ ! -s $simdir/$2/scores/all-alleles-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores.csv ]; then fail "Could not combine scores for iteration $i"; fi
	fi
	if [[ ! -s $simdir/$2/pr/all-alleles-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores-pr.csv || ! -s $simdir/$2/pr/all-alleles-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores-pr-each-labeled.csv ]]
	then
		$codedir/precision-recall.py -s $simdir/$2/scores/all-alleles-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores.csv | awk -v c="$2" -v e="$4" -v k="$7" -v i="$i" -v l="$3" -v f="$6" -v m="$8" -v OFS="," '{print $0, c, e, k, i, l, f, m}' > $simdir/$2/pr/all-alleles-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores-pr.csv
		$codedir/precision-recall-each-allele.py -s $simdir/$2/scores/all-alleles-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores.csv -o $simdir/$2/pr/all-alleles-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores-pr-each
		awk -v c="$2" -v e="$4" -v k="$7" -v i="$i" -v l="$3" -v f="$6" -v m="$8" -v OFS="," '{print $0, c, e, k, i, l, f, m}' $simdir/$2/pr/all-alleles-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores-pr-each.csv > $simdir/$2/pr/all-alleles-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores-pr-each-labeled.csv
		if [[ ! -s $simdir/$2/pr/all-alleles-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores-pr.csv || ! -s $simdir/$2/pr/all-alleles-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i-F-$7-k-$8-scores-pr-each-labeled.csv ]]; then fail "Scores could not be combined for iteration $i"; fi
	fi
	echo "Precision and recall for iteration $i complete"
done


# combine meta scores and determine how many times alleles called correctly across iterations
if [[ ! -s $simdir/pr-meta/all-alleles-all-i-flank10k-$3-bp-$6-frag-$2-x-e-$4-F-$7-k-$8-scores-pr.csv || ! -s $simdir/pr-meta/all-alleles-all-i-flank10k-$3-bp-$6-frag-$2-x-e-$4-F-$7-k-$8-scores-pr-each.csv ]]
then
	cat $simdir/$2/pr/all-alleles-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-*-F-$7-k-$8-scores-pr.csv > $simdir/pr-meta/all-alleles-all-i-flank10k-$3-bp-$6-frag-$2-x-e-$4-F-$7-k-$8-scores-pr.csv
	cat $simdir/$2/pr/all-alleles-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-*-F-$7-k-$8-scores-pr-each-labeled.csv > $simdir/pr-meta/all-alleles-all-i-flank10k-$3-bp-$6-frag-$2-x-e-$4-F-$7-k-$8-scores-pr-each.csv
	if [[ ! -s $simdir/pr-meta/all-alleles-all-i-flank10k-$3-bp-$6-frag-$2-x-e-$4-F-$7-k-$8-scores-pr.csv || ! -s $simdir/pr-meta/all-alleles-all-i-flank10k-$3-bp-$6-frag-$2-x-e-$4-F-$7-k-$8-scores-pr-each.csv ]]; then fail "Could not combine meta scores"; fi
fi

echo "Success!"
