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

# simulate reads noflank
wgsim -e 0 -d 250 -s 50 -N $num -1 100 -2 100 -r 0 -R 0 -X 0 -S $seed $alleledir/$a.fa $readdir/$a-noflank-100-bp-250-frag-20-x-e-0-i-$i.r1.reads.fq $readdir/$a-noflank-100-bp-250-frag-20-x-e-0-i-$i.r2.reads.fq
if [[ ! -s $readdir/$a-noflank-100-bp-250-frag-20-x-e-0-i-$i.r1.reads.fq || ! -s $readdir/$a-noflank-100-bp-250-frag-20-x-e-0-i-$i.r2.reads.fq ]]; then fail "Reads not simulated for noflank allele $a iteration $i"; fi

# cat reads
cat $readdir/$a-noflank-100-bp-250-frag-20-x-e-0-i-$i.r1.fq $readdir/$a-noflank-100-bp-250-frag-20-x-e-0-i-$i.r2.fq | paste - - - - | cut -f2 > $readdir/$a-noflank-100-bp-250-frag-20-x-e-0-i-$i.r12.reads

# run gbkc
$gbkcdir/gbkc count -a $thesisdir/PRDM9-36/all-alleles.fa -f $thesisdir/PRDM9-36/all-alleles.fa -1 $readdir/$a-noflank-100-bp-250-frag-20-x-e-0-i-$i.r1.reads.fq -2 $readdir/$a-noflank-100-bp-250-frag-20-x-e-0-i-$i.r2.reads.fq -K $k -l 100 -e 0 -c 20 -m coverage -o $savedir/gbkc/$a-noflank-100-bp-250-frag-20-x-e-0-i-$i-reads.scores
if [[ ! -s $savedir/gbkc/$a-noflank-100-bp-250-frag-20-x-e-0-i-$i-reads.scores ]]; then fail "gbkc failed for noflank allele $a iteration $i"; fi

echo "Reads simulated for iter $i"
done

echo success!

