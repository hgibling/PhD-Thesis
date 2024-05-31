#!/bin/bash
#$ -cwd

# command use: ./gbkc-reads.sh $1 allele list $2 coverage $3 read length $4 error rate $5 number of iterations $6 fragment length

# load modules
module load wgsim
module load python
source ~/scripts/error-handling.sh

# define directories
basedir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM
codedir=/.mounts/labs/simpsonlab/users/hgibling/BinfTK
gbkcdir=$basedir/gbkc/src
alleledir=$basedir/indiv-alleles
alleleflankdir=$alleledir/flank10k

simdir=$basedir/simulations/READS/haploid

# make directories if they don't already exist
mkdir -p $simdir/$2/reads


# start seed
seed=0


# iterate over list of alleles
a=$1
	seed=$((seed + 1))
	# calculate number of reads needed to be simulated for coverage
	num="$($codedir/get-read-number.py -s $alleleflankdir/$a-flank10k.fa.seq -c $2 -r $3 -p)"
	echo "$num reads"
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


echo "Success!"
