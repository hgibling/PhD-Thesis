#!/bin/bash
#$ -cwd
#$ -t 1-51

source ~/scripts/error-handling.sh
source /.mounts/labs/simpsonlab/users/hgibling/venvs/gbkc/bin/activate
module load python

basedir=/.mounts/labs/simpsonlab/users/hgibling/GQ-OHS
prdm9dir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM
gblrdir=/.mounts/labs/awadallalab/private/hgibling/gblr

sample=$(awk "NR==$SGE_TASK_ID" $basedir/OHS-samples.txt)

time $gblrdir/gblr.py -a $prdm9dir/all-alleles-standard-flank10k.fa -r $basedir/mapped-GRCh38-filtered/$sample-CCS-targeted-PRDM9-gap-filtered.bam -d -N 10 -c -e 0.01 -V -o $basedir/gblr-calls/$sample-gap-filtered-gblr.tsv
if [ $? -ne 0 ]; then fail "gblr failed"; fi

echo success!

