#!/bin/bash
#$ -cwd
#$ -t 1-49

basedir=/.mounts/labs/simpsonlab/users/hgibling/GQ-OHS
savedir=$basedir/gblr-calls-consensus-only
gblrdir=/.mounts/labs/awadallalab/private/hgibling/gblr
prdm9dir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM

source ~/scripts/error-handling.sh
source /.mounts/labs/simpsonlab/users/hgibling/venvs/gbkc/bin/activate
module load python

sample=$(awk "NR==$SGE_TASK_ID" $basedir/OHS-sample-names.txt)

$gblrdir/gblr.py -a $prdm9dir/all-alleles-standard-flank10k.fa -r $basedir/mapped-GRCh38-filtered-maxima/$sample-CCS-targeted-PRDM9-gap-filtered-maxima.bam -d -c -e 0.01 -s 1ee -E 1Indel -x consensus -o $savedir/$sample-gblr-GRCh38-consensus-only-f25.tsv
if [ $? -ne 0 ]; then fail "gblr failed"; fi

echo success!
