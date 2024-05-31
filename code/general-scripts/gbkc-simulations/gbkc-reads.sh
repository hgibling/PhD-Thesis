#!/bin/bash
#$ -cwd

# command use: ./gbkc-reads.sh $1 allele list $2 coverage $3 read length $4 error rate $5 number of iterations $6 fragment length $7 save directory

# load modules
module load wgsim
module load python/2.7.11
source /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/hmm-venv/bin/activate
source ~/scripts/error-handling.sh

# define directories
basedir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM
codedir=$basedir/PRDM9-Allele-Calling/code
gbkcdir=$basedir/gbkc/src
alleledir=$basedir/indiv-alleles
alleleflankdir=$alleledir/flank/flank10k

mkdir -p $basedir/simulations/$7
simdir=$basedir/simulations/$7


# make directories if they don't already exist
mkdir -p $simdir/$2/reads


# start seed
seed=0


# iterate over list of alleles
while read a
do
	seed=$((seed + 1))
	# calculate number of reads needed to be simulated for coverage
	num="$($codedir/calculate-coverage.py $alleleflankdir/$a-flank10k.csv $2 $3 T)"
	if [ -z $num ]; then fail "Number of reads to be simulated not calculated for allele $a"; fi

	for i in $(seq 1 $5)
	do
		if [ ! -s $simdir/$2/reads/$a-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i.fq ]; then
			seed=$((seed + 1))
			# simulate reads
			wgsim -e $4 -d $6 -s 50 -N $num -1 $3 -2 $3 -r 0 -R 0 -X 0 -S $seed $alleleflankdir/$a-flank10k.fa $simdir/$2/reads/$a-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i.r1.fq $simdir/$2/reads/$a-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i.r2.fq
			if [[ ! -s $simdir/$2/reads/$a-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i.r1.fq || ! -s $simdir/$2/reads/$a-flank10k-$3-bp-$6-frag-$2-x-e-$4-i-$i.r2.fq ]]; then fail "Reads not simulated for allele $a iteration $i"; fi
		fi
	done
	echo "Reads simulated for allele $a"
done < $1


echo "Success!"
