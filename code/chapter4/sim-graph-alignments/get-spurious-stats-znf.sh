#!/bin/bash
#$ -cwd
source ~/scripts/error-handling.sh

thesisdir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/THESIS
scriptdir=$thesisdir/PRDM9-36/graphs/scripts

a=$(awk "NR==$SGE_TASK_ID" $thesisdir/PRDM9-36/all-allele-names.txt)

module load rstats
Rscript --vanilla $scriptdir/spurious-alignments-znf.R $a

if [ $? -ne 0 ]; then fail "R failed"; fi
echo "success!"
