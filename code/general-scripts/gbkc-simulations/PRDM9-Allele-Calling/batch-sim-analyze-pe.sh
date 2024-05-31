#!/bin/bash
#$ -cwd

# command use: ./batch-sim-analyze.sh $1 allele list $2 coverage $3 read length $4 error rate $5 number of iterations $6 kmer size

# load modules
module load vg
module load python/2.7.11
source /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/hmm-venv/bin/activate
source ~/scripts/error-handling.sh


# define directories
basedir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM
codedir=$basedir/PRDM9-Allele-Calling/code
alleledir=$basedir/indiv-alleles
alleleflankdir=$alleledir/flank/flank1000
allelekmerdir=$basedir/kmers-distance
simdir=$basedir/sim-reads/pe-lam


# make directories if they don't already exist
mkdir -p $simdir/$2
mkdir -p $simdir/$2/reads
mkdir -p $simdir/$2/kmers
mkdir -p $simdir/$2/scores
mkdir -p $simdir/$2/analyzed
mkdir -p $simdir/$2/pr
mkdir -p $simdir/meta
mkdir -p $simdir/pr-meta


# generate kmer distances from alleles if they don't already exist
if [ ! -d $allelekmerdir/$6 ]
then
	mkdir $basedir/kmers-pos/$6 $allelekmerdir/$6
	while read a
	do
		$codedir/get-kmers.py $alleledir/$a.seq.csv $6 T > $basedir/kmers-pos/$6/$a-$6-k.kmerpos
		$codedir/kmer-distances.py -k $basedir/kmers-pos/$6/$a-$6-k.kmerpos -o $allelekmerdir/$6/$a-$6-k.kmerdistances
		if [ ! -s $allelekmerdir/$6/$a-$6-k.kmerdistances ]; then fail "No kmer distances generated for allele $a"; fi
		echo "Allele $a kmer distances complete"
		awk -v a="$a" -v OFS="," '{print a, $0}' $allelekmerdir/$6/$a-$6-k.kmerdistances >> $allelekmerdir/$6/all-alleles-$6-k.kmerdistances
	done < $1
fi

if [ ! -s $allelekmerdir/$6/all-alleles-$6-k.kmerdistances ]; then fail "No allele kmers found"; fi

# calculate lambda for scoring
# (readlength - k + 1) * (coverage / readlength)

# lam=$($codedir/calc-lambda.py $3 $6 $2 $4 1)
# if [ -z "$lam" ]; then fail "Lambda not calculated"; else echo "Using lambda value of $lam"; fi


# iterate over list of alleles
while read a
do
	# calculate number of reads needed to be simulated for coverage, if necessary
	num="$($codedir/coverage-calc.py $alleleflankdir/$a-flank1000.csv $2 $3)"
	if [ -z $num ]; then fail "Number of reads to be simulated not calculated for allele $a iteration $i"; fi

	for i in $(seq 1 $5)
	do
		# check if reads have already been simulated
		if [[ ! -s $simdir/$2/reads/$a-flank1000-$3-bp-$2-x-e-$4-i-$i-F-pe1.reads || ! -s $simdir/$2/reads/$a-flank1000-$3-bp-$2-x-e-$4-i-$i-F-pe2.reads ]]
		then
			# simulate reads
			vg sim -x $alleleflankdir/graphs/$a.xg -l $3 -n $num -f -e $4 -p 250 -v 50 > $simdir/$2/reads/$a-flank1000-$3-bp-$2-x-e-$4-i-$i-F-pe.reads
			cut -f1 $simdir/$2/reads/$a-flank1000-$3-bp-$2-x-e-$4-i-$i-F-pe.reads | awk -v a="$a" -v OFS="," '{print a "_1" ":" NR, $0}' > $simdir/$2/reads/$a-flank1000-$3-bp-$2-x-e-$4-i-$i-F-pe1.reads
			cut -f2 $simdir/$2/reads/$a-flank1000-$3-bp-$2-x-e-$4-i-$i-F-pe.reads | $codedir/reverse-complement.py - | awk -v a="$a" -v OFS="," '{print a "_2" ":" NR, $0}' > $simdir/$2/reads/$a-flank1000-$3-bp-$2-x-e-$4-i-$i-F-pe2.reads
			if [[ ! -s $simdir/$2/reads/$a-flank1000-$3-bp-$2-x-e-$4-i-$i-F-pe1.reads || ! -s $simdir/$2/reads/$a-flank1000-$3-bp-$2-x-e-$4-i-$i-F-pe2.reads ]]; then fail "Reads not simulated for allele $a iteration $i"; fi
			echo "Reads simulated for allele $a iteration $i complete"
		fi
		
		# check if kmers have already been generated for simulated reads
		if [[ ! -s $simdir/$2/kmers/$a-flank1000-$3-bp-$2-x-e-$4-i-$i-F-pe1-$6-k.kmerpos || ! -s $simdir/$2/kmers/$a-flank1000-$3-bp-$2-x-e-$4-i-$i-F-pe2-$6-k.kmerpos ]]
		then
			# get kmers
			$codedir/get-kmers-readID.py $simdir/$2/reads/$a-flank1000-$3-bp-$2-x-e-$4-i-$i-F-pe1.reads $6 T > $simdir/$2/kmers/$a-flank1000-$3-bp-$2-x-e-$4-i-$i-F-pe1-$6-k.kmerpos
			$codedir/get-kmers-readID.py $simdir/$2/reads/$a-flank1000-$3-bp-$2-x-e-$4-i-$i-F-pe2.reads $6 T > $simdir/$2/kmers/$a-flank1000-$3-bp-$2-x-e-$4-i-$i-F-pe2-$6-k.kmerpos
			if [[ ! -s $simdir/$2/kmers/$a-flank1000-$3-bp-$2-x-e-$4-i-$i-F-pe1-$6-k.kmerpos || ! -s $simdir/$2/kmers/$a-flank1000-$3-bp-$2-x-e-$4-i-$i-F-pe2-$6-k.kmerpos ]]; then fail "No kmers for allele $a iteration $i"; fi
			echo "kmers generated for allele $a iteration $i reads complete"
		fi
		
		# calculate scores
		$codedir/score-kmer-distances.py -a $allelekmerdir/$6/all-alleles-$6-k.kmerdistances -k1 $simdir/$2/kmers/$a-flank1000-$3-bp-$2-x-e-$4-i-$i-F-pe1-$6-k.kmerpos -k2 $simdir/$2/kmers/$a-flank1000-$3-bp-$2-x-e-$4-i-$i-F-pe2-$6-k.kmerpos > $simdir/$2/scores/$a-flank1000-$3-bp-$2-x-e-$4-i-$i-F-pe-$6-k-scores.csv
		if [ ! -s $simdir/$2/scores/$a-flank1000-$3-bp-$2-x-e-$4-i-$i-F-pe-$6-k-scores.csv ]; then fail "No scores for allele $a iteration $i"; fi
		echo "Scores for allele $a iteration $i complete"
	done
done < $1


# combine scores and analyze how many alleles called correctly for each set of conditions
for i in $(seq 1 $5)
do
	cat $simdir/$2/scores/*flank1000-$3-bp-$2-x-e-$4-i-$i-F-pe-$6-k-scores.csv > $simdir/$2/scores/all-alleles-flank1000-$3-bp-$2-x-e-$4-i-$i-F-pe-$6-k-scores.csv
	if [ ! -s $simdir/$2/scores/all-alleles-flank1000-$3-bp-$2-x-e-$4-i-$i-F-pe-$6-k-scores.csv ]; then fail "Could not combine scores for iteration $i"; fi

	$codedir/analyze-scores.py -s $simdir/$2/scores/all-alleles-flank1000-$3-bp-$2-x-e-$4-i-$i-F-pe-$6-k-scores.csv -o $simdir/$2/analyzed/all-alleles-flank1000-$3-bp-$2-x-e-$4-i-$i-F-pe-$6-k-scores-analyzed
	if [ ! -s $simdir/$2/analyzed/all-alleles-flank1000-$3-bp-$2-x-e-$4-i-$i-F-pe-$6-k-scores-analyzed.csv ]; then fail "Could not analyze scores for iteration $i"; fi
	echo "Correct number of alleles called for iteration $i complete"
done


# calculate precision and recall overall and for each allele
for i in $(seq 1 $5)
do
	$codedir/pr-analyze-scores.py -s $simdir/$2/scores/all-alleles-flank1000-$3-bp-$2-x-e-$4-i-$i-F-pe-$6-k-scores.csv | awk -v c="$2" -v e="$4" -v k="$6" -v i="$i" -v OFS="," '{print $0, c, e, k, i}' > $simdir/$2/pr/all-alleles-flank1000-$3-bp-$2-x-e-$4-i-$i-F-pe-$6-k-scores-pr.csv
	$codedir/pr-analyze-scores-each.py -s $simdir/$2/scores/all-alleles-flank1000-$3-bp-$2-x-e-$4-i-$i-F-pe-$6-k-scores.csv -o $simdir/pr/$simdir/$2/pr/all-alleles-flank1000-$3-bp-$2-x-e-$4-i-$i-F-pe-$6-k-scores-pr-each
	awk -v c="$2" -v e="$4" -v k="$6" -v i="$i" -v OFS="," '{print $0, c, e, k, i}' $simdir/pr/$simdir/$2/pr/all-alleles-flank1000-$3-bp-$2-x-e-$4-i-$i-F-pe-$6-k-scores-pr-each.csv > $simdir/$2/pr/all-alleles-flank1000-$3-bp-$2-x-e-$4-i-$i-F-pe-$6-k-scores-pr-each-labeled.csv
	if [[ ! -s $simdir/$2/pr/all-alleles-flank1000-$3-bp-$2-x-e-$4-i-$i-F-pe-$6-k-scores-pr.csv || ! -s $simdir/$2/pr/all-alleles-flank1000-$3-bp-$2-x-e-$4-i-$i-F-pe-$6-k-scores-pr.csv ]]; then fail "Scores could not be combined for iteration $i"; fi
	echo "Precision and recall for iteration $i complete"
done

# # combine meta scores and determine how many times alleles called correctly across iterations
# cat $simdir/$2/analyzed/all-alleles-flank1000-$3-bp-$2-x-e-$4-i-*-F-$6-k-scores-analyzed.csv | sort -n | uniq -c | awk -v i="$2" -v j="$4" -v k="$6" -v OFS="," '{print $1, i, j, k, $2}' > $simdir/meta/meta-flank1000-$3-bp-$2-x-e-$4-F-pe-$6-k.csv
# cat $simdir/$2/pr/all-alleles-flank1000-$3-bp-$2-x-e-$4-i-*-F-pe-$6-k-scores-pr.csv > $simdir/$2/pr-meta/all-alleles-all-i-flank1000-$3-bp-$2-x-e-$4-F-pe-$6-k-scores-pr.csv
# cat $simdir/$2/pr/all-alleles-flank1000-$3-bp-$2-x-e-$4-i-*-F-pe-$6-k-scores-pr-each-labeled.csv > $simdir/$2/pr-meta/all-alleles-all-i-flank1000-$3-bp-$2-x-e-$4-F-pe-$6-k-scores-pr-each.csv
# if [[ ! -s $simdir/meta/meta-flank1000-$3-bp-$2-x-e-$4-F-pe-$6-k.csv || ! -s $simdir/$2/pr-meta/all-alleles-all-i-flank1000-$3-bp-$2-x-e-$4-F-pe-$6-k-scores-pr.csv || ! -s $simdir/$2/pr-meta/all-alleles-all-i-flank1000-$3-bp-$2-x-e-$4-F-pe-$6-k-scores-pr-each.csv ]]; then fail "Could not combine meta scores"; fi


echo "Success!"