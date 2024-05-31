#!/bin/bash
#$ -cwd
#$ -t 1-666
g=$(awk "NR==$SGE_TASK_ID" /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/diploid-pairs.txt)
/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/scripts/gbkc-distances-diploid-part1.sh $g 40 100 0.001 100 250 99 max distances-diploid
