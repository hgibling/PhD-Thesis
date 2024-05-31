#!/bin/bash
#$ -cwd
#$ -t 1-53
source ~/scripts/error-handling.sh
sample=$(awk "NR==$SGE_TASK_ID" /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/validation/long-read-calling/sample-names.txt)
source /.mounts/labs/simpsonlab/users/hgibling/venvs/gbkc/bin/activate
module load python mafft
/.mounts/labs/awadallalab/private/hgibling/gblr/gblr.py -a /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/all-alleles-standard-flank10k.fa -r /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/validation/long-read-calling/GRCh38-bams/$sample-ccs-prdm9-10k-GRCh38.bam -d -N 10 -v -V -R chr5:23526673,23527764 -e 0.01 -o /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/validation/long-read-calling/gblr-calls/$sample-gblr-GRCh38.tsv
if [ $? -ne 0 ]; then fail "gblr failed"; fi
