#!/bin/bash
#$ -cwd
source ~/scripts/error-handling.sh
/.mounts/labs/awadallalab/private/hgibling/gbkc/src/gbkc check-profiles -a /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/all-alleles-flank10k.fa -l 201 -u 300 -d > /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/count-profiles/diploid-profiles-300k.csv
