#!/bin/bash
#$ -cwd
#$ -t 1-666
g=$(awk "NR==$SGE_TASK_ID" /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/diploid-pairs.txt)
/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/scripts/gbkc-reads-diploid-simpson-fake-200x-flank1k.sh $g 200 100 0 25 250 200x-actually-flank1k
