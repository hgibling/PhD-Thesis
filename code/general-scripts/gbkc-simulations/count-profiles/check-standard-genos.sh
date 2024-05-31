#!/bin/bash
#$ -cwd
#$ -t 1-10
source ~/scripts/error-handling.sh

gbkcdir=/.mounts/labs/awadallalab/private/hgibling/gbkc/src
prdm9dir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM
savedir=$prdm9dir/count-profiles

$gbkcdir/gbkc check-profiles -k 300 -K 310 -d -a $prdm9dir/all-alleles-standard-flank10k.fa -o $savedir/standard-allele-genotype-unique.txt -t 20

echo "success!"
