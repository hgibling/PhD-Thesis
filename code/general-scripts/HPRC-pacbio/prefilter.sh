#!/bin/bash
#$ -cwd
#$ -t 1-52

basedir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/validation/long-read-calling
savedir=$basedir/prefiltered-bams
gblrdir=/.mounts/labs/awadallalab/private/hgibling/gblr
prdm9dir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM

source ~/scripts/error-handling.sh
source /.mounts/labs/simpsonlab/users/hgibling/venvs/gbkc/bin/activate
module load python

sample=$(awk "NR==$SGE_TASK_ID" $basedir/sample-names-used.txt)

$gblrdir/gap-filter.py -b $basedir/GRCh38-bams/$sample-ccs-prdm9-10k-GRCh38.bam -g 100000 -R 0 -v -o $savedir/$sample-keep-reads.txt
if [ $? -ne 0 ]; then fail "gblr failed"; fi

echo success!
