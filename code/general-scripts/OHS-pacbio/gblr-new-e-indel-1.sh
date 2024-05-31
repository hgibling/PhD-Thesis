#!/bin/bash
#$ -cwd
#$ -t 1-50

basedir=/.mounts/labs/simpsonlab/users/hgibling/GQ-OHS/
bamdir=$basedir/mapped-GRCh38-filtered-maxima
savedir=$basedir/gblr-calls-new-e-indel-1
gblrdir=/.mounts/labs/awadallalab/private/hgibling/gblr
prdm9dir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM

sample=$(awk "NR==$SGE_TASK_ID" $basedir/OHS-sample-names.txt)

source ~/scripts/error-handling.sh
source /.mounts/labs/simpsonlab/users/hgibling/venvs/gbkc/bin/activate
module load python

$gblrdir/gblr.py -a $prdm9dir/all-alleles-standard-flank10k.fa -r $bamdir/$sample-CCS-targeted-PRDM9-gap-filtered-maxima.bam -d -N 10 -c -e 0.01 -s 1ee -E 1Indel -o $savedir/$sample-gblr-GRCh38.tsv
if [ $? -ne 0 ]; then fail "gblr failed"; fi

echo success!
