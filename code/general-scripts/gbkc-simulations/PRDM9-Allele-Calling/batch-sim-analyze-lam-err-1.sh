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
allelekmerdir=$basedir/kmers
simdir=$basedir/sim-reads/lam-err-redo


# make directories if they don't already exist
mkdir -p $simdir/$2
mkdir -p $simdir/$2/reads
mkdir -p $simdir/$2/kmers
mkdir -p $simdir/$2/scores
mkdir -p $simdir/$2/analyzed
mkdir -p $simdir/$2/pr
mkdir -p $simdir/meta
mkdir -p $simdir/pr-meta


# # generate kmers from alleles if they don't already exist
# if [ ! -d $allelekmerdir/$6 ]
# then
# 	mkdir $allelekmerdir/$6
# 	while read a
# 	do
# 		$codedir/get-kmers.py $alleledir/$a.seq.csv $6 | sort | uniq -c | awk -v a="$a" '{print a "," $2 "," $1}' > $allelekmerdir/$6/$a-$6-k.kmercounts
# 		if [ ! -s $allelekmerdir/$6/$a-$6-k.kmercounts ]; then fail "No kmers generated for allele $a"; fi
# 	done < $1
# 	cat $allelekmerdir/$6/*-$6-k.kmercounts > $allelekmerdir/$6/all-alleles-$6-k.kmercounts
# fi

# if [ ! -s $allelekmerdir/$6/all-alleles-$6-k.kmercounts ]; then fail "No allele kmers found"; fi


# calculate lambda for scoring
# (readlength - k + 1) * (coverage / readlength)

lam=$($codedir/calc-lambda.py $3 $6 $2 $4 1)
if [ -z "$lam" ]; then fail "Lambda not calculated"; else echo "Using lambda value of $lam"; fi


# iterate over list of alleles
while read a
do
	for i in $(seq 1 $5)
	do
		# # check if reads have already been simulated
		# if [ ! -f $simdir/$2/reads/$a-flank1000-$3-bp-$2-x-e-$4-i-$i-F.reads ]
		# then
		# 	# calculate number of reads needed to be simulated for coverage
		#         num="$($codedir/coverage-calc.py $alleleflankdir/$a-flank1000.csv $2 $3)"
		# 	if [ -z $num ]; then fail "Number of reads to be simulated not calculated for allele $a iteration $i"; fi

		# 	# simulate reads
		# 	vg sim -x $alleleflankdir/graphs/$a.xg -l $3 -n $num -f -e $4 > $simdir/$2/reads/$a-flank1000-$3-bp-$2-x-e-$4-i-$i-F.reads
		# 	if [ ! -s $simdir/$2/reads/$a-flank1000-$3-bp-$2-x-e-$4-i-$i-F.reads ]; then fail "Reads not simulated for allele $a iteration $i"; fi
		# fi

		# check if kmers have already been generated for simulated reads
		if [ ! -s $simdir/$2/kmers/$a-flank1000-$3-bp-$2-x-e-$4-i-$i-F-$6-k.kmercounts ]
		then
			# get kmers
			$codedir/get-kmers.py $simdir/$2/reads/$a-flank1000-$3-bp-$2-x-e-$4-i-$i-F.reads $6 | sort | uniq -c | awk -v a="$a" '{print a "," $2 "," $1}' > $simdir/$2/kmers/$a-flank1000-$3-bp-$2-x-e-$4-i-$i-F-$6-k.kmercounts
			if [ ! -s $simdir/$2/kmers/$a-flank1000-$3-bp-$2-x-e-$4-i-$i-F-$6-k.kmercounts ]; then fail "No kmers for allele $a iteration $i"; fi
		fi

		# calculate scores
		$codedir/tb-count-score-all-alleles.py -r $simdir/$2/kmers/$a-flank1000-$3-bp-$2-x-e-$4-i-$i-F-$6-k.kmercounts -a $allelekmerdir/$6/all-alleles-$6-k.kmercounts -l $lam -e 1 > $simdir/$2/scores/$a-flank1000-$3-bp-$2-x-e-$4-i-$i-F-$6-k-scores.csv
		if [ ! -s $simdir/$2/scores/$a-flank1000-$3-bp-$2-x-e-$4-i-$i-F-$6-k-scores.csv ]; then fail "No scores for allele $a iteration $i"; fi
	done
done < $1


# combine scores and analyze how many alleles called correctly for each set of conditions
for i in $(seq 1 $5)
do
	cat $simdir/$2/scores/*flank1000-$3-bp-$2-x-e-$4-i-$i-F-$6-k-scores.csv > $simdir/$2/scores/all-alleles-flank1000-$3-bp-$2-x-e-$4-i-$i-F-$6-k-scores.csv
	if [ ! -s $simdir/$2/scores/all-alleles-flank1000-$3-bp-$2-x-e-$4-i-$i-F-$6-k-scores.csv ]; then fail "Could not combine scores for iteration $i"; fi

	$codedir/analyze-scores.py -s $simdir/$2/scores/all-alleles-flank1000-$3-bp-$2-x-e-$4-i-$i-F-$6-k-scores.csv -o $simdir/$2/analyzed/all-alleles-flank1000-$3-bp-$2-x-e-$4-i-$i-F-$6-k-scores-analyzed
	if [ ! -s $simdir/$2/analyzed/all-alleles-flank1000-$3-bp-$2-x-e-$4-i-$i-F-$6-k-scores-analyzed.csv ]; then fail "Could not analyze scores for iteration $i"; fi
done


# combine meta scores and determine how many times alleles called correctly across iterations
cat $simdir/$2/analyzed/all-alleles-flank1000-$3-bp-$2-x-e-$4-i-*-F-$6-k-scores-analyzed.csv | sort -n | uniq -c | awk -v i="$2" -v j="$4" -v k="$6" -v OFS="," '{print $1, i, j, k, $2}' > $simdir/meta/meta-flank1000-$3-bp-$2-x-e-$4-$6-k.csv
if [ ! -s $simdir/meta/meta-flank1000-$3-bp-$2-x-e-$4-$6-k.csv ]; then fail "Could not combine meta scores"; fi


echo "Success!"
