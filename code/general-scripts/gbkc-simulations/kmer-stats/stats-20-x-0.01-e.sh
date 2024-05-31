#!/bin/bash
#$ -cwd
#$ -t 1-23
k=$(awk "NR==$SGE_TASK_ID" /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/k-fewer.txt)
/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/scripts/gbkc-counts-kmer-stats.sh /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/names.txt 20 100 0.01 10 250 $k kmer-stats
