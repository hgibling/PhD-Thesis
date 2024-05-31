#!/bin/bash
#$ -cwd
#$ -t 1-52

basedir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/validation/long-read-calling
savedir=$basedir/gblr-calls-compare-consensus-seqs
gblrdir=/.mounts/labs/awadallalab/private/hgibling/gblr
prdm9dir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM

source ~/scripts/error-handling.sh
source /.mounts/labs/simpsonlab/users/hgibling/venvs/gbkc/bin/activate
module load python

sample=$(awk "NR==$SGE_TASK_ID" $basedir/sample-names-used.txt)

$gblrdir/gblr.py -a $prdm9dir/all-alleles-standard-flank10k.fa -r $basedir/GRCh38-bams/$sample-ccs-prdm9-10k-GRCh38.bam -d -N 10 -c -s 1e -E 1Indel -e 0.01 -o $savedir/$sample-gblr-GRCh38-1e-1i.tsv
if [ $? -ne 0 ]; then fail "gblr failed"; fi

echo success!
