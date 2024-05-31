#!/bin/bash
#$ -cwd
#$ -t 1-140
l=$(awk "NR==$SGE_TASK_ID" /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/simulations/manual-lambda/lambda-values.txt)
/.mounts/labs/awadallalab/private/hgibling/gbkc/src/gbkc count -a /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/all-alleles.fa -f /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/flank10k.fa -1 /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/simulations/READS/diploid/100/reads/A-A-flank10k-100-bp-250-frag-100-x-e-0-i-5.r1.fq -2 /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/simulations/READS/diploid/100/reads/A-A-flank10k-100-bp-250-frag-100-x-e-0-i-5.r2.fq -c 100 -e 0 -l 100 -k 51 -K 71 -i 20 -d -m coverage -M $l -o /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/simulations/manual-lambda/A-A-flank10k-100-bp-250-frag-100-x-e-0-i-5-lam-$l.csv
