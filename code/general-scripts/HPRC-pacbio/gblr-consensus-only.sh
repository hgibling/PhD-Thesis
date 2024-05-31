#!/bin/bash
#$ -cwd
#$ -t 1-52

basedir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/validation/long-read-calling
savedir=$basedir/gblr-calls-consensus-only
gblrdir=/.mounts/labs/awadallalab/private/hgibling/gblr
prdm9dir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM

source ~/scripts/error-handling.sh
source /.mounts/labs/simpsonlab/users/hgibling/venvs/gbkc/bin/activate
module load python

sample=$(awk "NR==$SGE_TASK_ID" $basedir/sample-names-used.txt)

$gblrdir/gblr.py -a $prdm9dir/all-alleles-standard-flank10k.fa -r $basedir/GRCh38-bams/$sample-ccs-prdm9-10k-GRCh38.bam -d  -c -e 0.01 -s 1ee -E 1Indel -x consensus -o $savedir/$sample-gblr-GRCh38-consensus-only.tsv
if [ $? -ne 0 ]; then fail "gblr failed"; fi

echo success!
