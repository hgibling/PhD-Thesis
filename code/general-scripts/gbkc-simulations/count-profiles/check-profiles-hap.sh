#!/bin/bash
#$ -cwd
source ~/scripts/error-handling.sh
/.mounts/labs/awadallalab/private/hgibling/gbkc/src/gbkc check-profiles -a /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/all-alleles-flank10k.fa -k 1 -K 100 > /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/count-profiles/hap-profiles-100k.csv
