#!/bin/bash
#$ -cwd
#$ -t 1-666
g=$(awk "NR==$SGE_TASK_ID" /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/diploid-pairs.txt)
/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/scripts/gbkc-counts-diploid-simpson-scratch-real-200x-part1-flank1k.sh $g 200 100 0 25 250 99 NA
