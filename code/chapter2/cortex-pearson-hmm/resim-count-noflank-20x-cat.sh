#!/bin/bash
#$ -cwd

# load modules
module load wgsim
module load python/2.7
source ~/scripts/error-handling.sh

basedir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM
codedir=$basedir/PRDM9-Allele-Calling/code
gbkcdir=/.mounts/labs/awadallalab/private/hgibling/gbkc/src
alleledir=$basedir/indiv-alleles
alleleflankdir=$alleledir/flank10k
thesisdir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/THESIS
readdir=/scratch2/groups/simpsonlab/hgibling/no-flank-20x
savedir=$thesisdir/PRDM9-36/HMM-test/no-flank-20x

a=$(awk "NR==$SGE_TASK_ID" $basedir/names.txt)
seed=0
k=99

# calculate number of reads needed to be simulated for coverage
num="$($codedir/calculate-coverage.py $alleledir/$a.csv 20 100 T)"
if [[ -z $num ]]; then fail "Number of reads to be simulated not calculated for allele $a"; fi

# iterations
for i in $(seq 1 50)
do
seed=$((seed + 1))

# cat reads
cat $readdir/$a-noflank-100-bp-250-frag-20-x-e-0-i-$i.r1.reads.fq $readdir/$a-noflank-100-bp-250-frag-20-x-e-0-i-$i.r2.reads.fq | paste - - - - | cut -f2 > $readdir/$a-noflank-100-bp-250-frag-20-x-e-0-i-$i.r12.reads

echo "Reads simulated for iter $i"
done

echo success!

