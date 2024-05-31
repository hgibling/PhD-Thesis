#!/bin/bash
#$ -cwd

# command use: ./gbkc-reads.sh $1 genotype $2 coverage $3 read length $4 error rate $5 number of iterations $6 fragment length $7 save directory

# load modules
module load wgsim
module load python/2.7.11
module load jq
source /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/hmm-venv/bin/activate
source ~/scripts/error-handling.sh

# define directories
basedir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM
codedir=$basedir/PRDM9-Allele-Calling/code
gbkcdir=$basedir/gbkc/src
genotypedir=$basedir/indiv-genotypes/flank10k
alleledir=$basedir/indiv-alleles/flank10k

mkdir -p /.mounts/labs/simpsonlab/users/hgibling/simulations/$7
simdir=/.mounts/labs/simpsonlab/users/hgibling/simulations/$7


# make directories if they don't already exist
mkdir -p $simdir/100/reads


# start seed
seed=0


# iterate over list of genotypes
g=$1
	seed=$((seed + 1))

	# get individual alleles
	a1=$(echo $g | cut -d- -f1)
    a2=$(echo $g | cut -d- -f2)

	# calculate number of reads needed to be simulated for coverage
	num1="$($codedir/calculate-coverage.py $alleledir/$a1-flank10k.fa.csv $2 $3 T)"
	num2="$($codedir/calculate-coverage.py $alleledir/$a2-flank10k.fa.csv $2 $3 T)"
	sum=$((num1 + num2))
	num=$(jq -n $sum/2)

	if [ -z $num ]; then fail "Number of reads to be simulated not calculated for genotype $g"; fi

	for i in $(seq 1 $5)
	do
		if [[ ! -s $simdir/100/reads/$g-flank10k-$3-bp-$6-frag-100-x-e-$4-i-$i.r1.fq || ! -s $simdir/100/reads/$g-flank10k-$3-bp-$6-frag-100-x-e-$4-i-$i.r2.fq ]]; then
			seed=$((seed + 1))
			# simulate reads
			wgsim -e $4 -d $6 -s 50 -N $num -1 $3 -2 $3 -r 0 -R 0 -X 0 -S $seed $genotypedir/$g-flank10k.fa $simdir/100/reads/$g-flank10k-$3-bp-$6-frag-100-x-e-$4-i-$i.r1.fq $simdir/100/reads/$g-flank10k-$3-bp-$6-frag-100-x-e-$4-i-$i.r2.fq
			if [[ ! -s $simdir/100/reads/$g-flank10k-$3-bp-$6-frag-100-x-e-$4-i-$i.r1.fq || ! -s $simdir/100/reads/$g-flank10k-$3-bp-$6-frag-100-x-e-$4-i-$i.r2.fq ]]; then fail "Reads not simulated for genotype $g iteration $i"; fi
		fi
		echo "Reads simulated for iteration $i"

	done



echo "Success!"
