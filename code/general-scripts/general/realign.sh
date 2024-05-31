#!/bin/bash
#$ -cwd
source ~/scripts/error-handling.sh

module load bwa
savedir=/.mounts/labs/simpsonlab/users/hgibling/spades-test/realign
refdir=/.mounts/labs/simpsonlab/users/hgibling/references

bwa mem -t 4 $refdir/hs37d5.fa $savedir/HG003.2x250_1.paired.fq $savedir/HG003.2x250_2.paired.fq > $savedir/HG003.2x250_original.sam

