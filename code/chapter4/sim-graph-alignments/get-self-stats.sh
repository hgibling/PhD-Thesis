#!/bin/bash
#$ -cwd
#$ -t 1-36
a=$(awk "NR==$SGE_TASK_ID" /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/names.txt)
module load rstats
Rscript --vanilla /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/graphs/scripts/analysis/pull-reads.R $a
