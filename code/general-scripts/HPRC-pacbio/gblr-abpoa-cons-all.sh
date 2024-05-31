#!/bin/bash
#$ -cwd
#$ -t 1-53

basedir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/validation/long-read-calling
prdm9dir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM
gblrdir=/.mounts/labs/awadallalab/private/hgibling/gblr

source ~/scripts/error-handling.sh
source /.mounts/labs/simpsonlab/users/hgibling/venvs/gbkc/bin/activate
sample=$(awk "NR==$SGE_TASK_ID" $basedir/sample-names-used.txt)
module load python

$gblrdir/gblr.py -a $prdm9dir/all-alleles-standard-flank10k.fa -r $basedir/GRCh38-bams/$sample-ccs-prdm9-10k-GRCh38.bam -d -o $basedir/gblr-calls-abpoa-all/$sample-gblr-GRCh38-abpoa.tsv
if [ $? -ne 0 ]; then fail "gblr failed"; fi

echo "success!"
