#!/bin/bash
#$ -cwd
#$ -t 1-52
source ~/scripts/error-handling.sh

module load rstats/4.2.2

paperdir=/.mounts/labs/simpsonlab/users/hgibling/PAPER
savedir=$paperdir/OHS/illumina
scriptdir=$paperdir/scripts

sample=$(awk "NR==$SGE_TASK_ID" $savedir/OHS-sample-names.txt | cut -f1)

Rscript --vanilla OHS-combine-scores-uncorrected-106.R $sample

if [ $? -ne 0 ]; then fail "combine failed"; fi
echo "success!"
