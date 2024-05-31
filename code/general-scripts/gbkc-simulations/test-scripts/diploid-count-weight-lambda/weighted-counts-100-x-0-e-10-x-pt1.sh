#!/bin/bash
#$ -cwd
#$ -t 1-666
g=$(awk "NR==$SGE_TASK_ID" /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/diploid-pairs.txt)
/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/scripts/gbkc-counts-diploid-simpson-scratch-weighted-part1.sh $g 100 100 0 25 250 99 10
