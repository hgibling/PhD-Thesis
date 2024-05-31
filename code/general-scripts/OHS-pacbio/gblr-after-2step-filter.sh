#!/bin/bash
#$ -cwd
#$ -t 1-50

basedir=/.mounts/labs/simpsonlab/users/hgibling/GQ-OHS/
bamdir=$basedir/mapped-GRCh38-filtered-maxima
savedir=$basedir/gblr-calls-2step-filtered
gblrdir=/.mounts/labs/awadallalab/private/hgibling/gblr
prdm9dir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM

source ~/scripts/error-handling.sh
sample=$(awk "NR==$SGE_TASK_ID" $basedir/OHS-samples.txt)
source /.mounts/labs/simpsonlab/users/hgibling/venvs/gbkc/bin/activate
module load python

$gblrdir/gblr.py -a $prdm9dir/all-alleles-standard-flank10k.fa -r $bamdir/$sample-CCS-targeted-PRDM9-gap-filtered-maxima.bam -d -N 10 -c -e 0.01 -V -o $savedir/$sample-2step-filtered-gblr-GRCh38.tsv
if [ $? -ne 0 ]; then fail "gblr failed"; fi

echo "Success!"
