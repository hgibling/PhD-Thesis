#!/bin/bash
#$ -cwd

module load python/2.7

source ~/scripts/error-handling.sh
source /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/thesispy2/bin/activate

cd /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/THESIS/PRDM9-36/HMM-test

python call-alleles.py -a A-short.alleles  -r A-short.reads -e 0.0000000001 -o A-test-rc-print2
if [ $? -ne 0 ]; then fail "fail"; fi
echo success!

